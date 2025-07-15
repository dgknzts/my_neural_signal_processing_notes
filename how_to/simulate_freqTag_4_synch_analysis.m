%% Simulating EEG data with phase-lagged dipoles
% Load empty EEG structure with leadfield
addpath("src");
load emptyEEG.mat

% Select two dipole locations
diploc1 = 101;
diploc2 = 201;

% EEG parameters
EEG.times = (-1:1/EEG.srate:2.4);
EEG.pnts = numel(EEG.times);
EEG.trials = 80;
EEG.data = zeros(EEG.nbchan, EEG.pnts, EEG.trials);

% Signal parameters
freq = 15; % Hz
phaselag = 0.15 * 2* pi; % phase difference between dipoles
time_range = [0 2]; % activation window

% Create rectangular window for signal activation
rect_win = (EEG.times >= time_range(1)) & (EEG.times <= time_range(2));

% Add some smoothing to avoid sharp edges
smooth_win = rect_win;
smooth_len = round(0.1 * EEG.srate); % 100ms smoothing
smooth_win = conv(double(smooth_win), hann(smooth_len)/sum(hann(smooth_len)), 'same');

% Initialize dipole activity matrix
dipole_data = zeros(size(lf.Gain,3), EEG.pnts);

%% Loop over trials
rng(1)

for triali = 1:EEG.trials
    
    % Reset dipole activity with background noise
    dipole_data = randn(size(dipole_data));
    
    % Add trial-specific phase jitter (optional - makes it more realistic)
    % trial_phase_jitter1 = 0.1 * randn; % small random phase variation
    % trial_phase_jitter2 = 0.1 * randn;
    
    % Create signals for both dipoles
    sig1 = sin(2*pi*freq*EEG.times) .* smooth_win;
    sig2 = sin(2*pi*freq*EEG.times + phaselag) .* smooth_win;
    
    % Add amaplitude variation across trials
    amp1 = 6 + 0.05*randn;
    amp2 = 6 + 0.05*randn;
    
    % Assign signals to dipoles
    dipole_data(diploc1,:) = amp1 * sig1;
    dipole_data(diploc2,:) = amp2 * sig2;
    
    % Project dipole activity to scalp electrodes using leadfield
    EEG.data(:,:,triali) = squeeze(lf.Gain(:,1,:)) * dipole_data;
end

% Add some pink noise to make it more realistic
for chani = 1:EEG.nbchan
    for triali = 1:EEG.trials
        % Generate pink noise
        ed = 50; % exponential decay
        as = rand(1,floor(EEG.pnts/2)-1) .* exp(-(1:floor(EEG.pnts/2)-1)/ed);
        as = [as(1) as 0 as(:,end:-1:1)];
        fc = as .* exp(1i*2*pi*rand(size(as)));
        noise = real(ifft(fc)) * EEG.pnts * 750; % scale noise
        
        % Add noise to signal
        EEG.data(chani,:,triali) = EEG.data(chani,:,triali) + noise;
    end
end
%% Show topographical maps of dipole projections
figure(1), clf
subplot(121)
topoplotIndie(-lf.Gain(:,1,diploc1), EEG.chanlocs, 'numcontour', 0);
title('Dipole 1 projection')
colorbar
subplot(122)
topoplotIndie(-lf.Gain(:,1,diploc2), EEG.chanlocs, 'numcontour', 0);
title('Dipole 2 projection')
colorbar

%% Check FFT spectrum for channels capturing each dipole
% Find channels with maximum projection from each dipole
[~, chan1] = max(abs(lf.Gain(:,1,diploc1)));
[~, chan2] = max(abs(lf.Gain(:,1,diploc2)));

% Extract data from activation window (0-2 seconds)
time_idx = (EEG.times >= 0) & (EEG.times <= 2);
n_points = sum(time_idx);

% Compute FFT for both channels
fft_chan1 = fft(mean(EEG.data(chan1, time_idx, :), 3), n_points);
fft_chan2 = fft(mean(EEG.data(chan2, time_idx, :), 3), n_points);

% Frequency vector
frequencies = linspace(0, EEG.srate, n_points);
freq_idx = frequencies <= 50; % only plot up to 50 Hz

% Compute power spectrum
pow_chan1 = 2*abs(fft_chan1(freq_idx))/n_points;
pow_chan2 = 2*abs(fft_chan2(freq_idx))/n_points;

figure(2), clf
subplot(211)
plot(frequencies(freq_idx), pow_chan1, 'b', 'LineWidth', 2)
hold on
plot([15 15], [0 max(pow_chan1)*1.1], 'r--', 'LineWidth', 2)
xlabel('Frequency (Hz)'), ylabel('Amplitude')
title(['FFT Spectrum - Channel ' num2str(chan1) ' (Best for Dipole 1)'])
xlim([0 50])
grid on

subplot(212)
plot(frequencies(freq_idx), pow_chan2, 'r', 'LineWidth', 2)
hold on
plot([15 15], [0 max(pow_chan2)*1.1], 'r--', 'LineWidth', 2)
xlabel('Frequency (Hz)'), ylabel('Amplitude')
title(['FFT Spectrum - Channel ' num2str(chan2) ' (Best for Dipole 2)'])
xlim([0 50])
grid on

%% Topoplot at 15 Hz with highest amplitude channel marked

% Compute FFT for all channels (average across trials, activation window only)
all_fft = zeros(EEG.nbchan, n_points);
for ch = 1:EEG.nbchan
    all_fft(ch,:) = abs(fft(mean(EEG.data(ch, time_idx, :), 3), n_points));
end

% Find the index closest to 15 Hz
[~, freq15_idx] = min(abs(frequencies - 15));

% Get amplitude at 15 Hz for all channels
amp_at_15Hz = 2*all_fft(:, freq15_idx)/n_points;

% Assuming amp_at_15Hz is a row vector
[sorted_amps, sorted_indices] = sort(amp_at_15Hz, 'descend');

figure(3), clf
topoplotIndie(amp_at_15Hz, EEG.chanlocs, 'numcontour', 0);
hold on
% Mark the highest amplitude channel
topoplotIndie(amp_at_15Hz, EEG.chanlocs, 'electrodes', 'on');
title(sprintf('Topography at 15 Hz\nMax channels: %d and %d', sorted_indices(1), sorted_indices(2)))
colorbar

save('G:\My Drive\Projects\signal_processing_mike_cohen/how_to/src/simulated_phaselag_EEG.mat', 'EEG', 'diploc1', 'diploc2', 'freq', 'phaselag');