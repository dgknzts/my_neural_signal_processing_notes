function EEG = simulate_phaseLag_data(freq, phaselag_rad, diploc1, diploc2, time_range, show_plots)
% simulate_eeg_data generates simulated EEG data from two phase-lagged dipoles.
%
% Inputs:
%   freq          - Signal frequency in Hz (e.g., 15)
%   phaselag_rad  - Phase difference between dipoles in radians (e.g., 0.125 * 2 * pi)
%   diploc1       - Index of the first dipole location (e.g., 101)
%   diploc2       - Index of the second dipole location (e.g., 201)
%   time_range    - 1x2 vector for signal activation window [start, end] in seconds (e.g., [0 2])
%   show_plots    - Boolean (true/false) to control plotting
%
% Output:
%   EEG           - EEG data structure containing the simulated data

%% Setup
% Load empty EEG structure with leadfield
addpath("src");
load emptyEEG.mat

% EEG parameters
EEG.times = (-1:1/EEG.srate:2.4);
EEG.pnts = numel(EEG.times);
EEG.trials = 80;
EEG.data = zeros(EEG.nbchan, EEG.pnts, EEG.trials);

% Create rectangular window for signal activation
rect_win = (EEG.times >= time_range(1)) & (EEG.times <= time_range(2));

% Add some smoothing to avoid sharp edges
smooth_len = round(0.1 * EEG.srate); % 100ms smoothing
smooth_win = conv(double(rect_win), hann(smooth_len)/sum(hann(smooth_len)), 'same');

% Initialize dipole activity matrix
dipole_data = zeros(size(lf.Gain,3), EEG.pnts);

%% Data Simulation
rng(1); % for reproducibility
for triali = 1:EEG.trials
    
    % Reset dipole activity with background noise
    dipole_data = randn(size(dipole_data));
    
    % Create signals for both dipoles
    sig1 = sin(2*pi*freq*EEG.times) .* smooth_win;
    sig2 = sin(2*pi*freq*EEG.times + phaselag_rad) .* smooth_win;
    
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
        ed = 50; % exponential decay
        as = rand(1,floor(EEG.pnts/2)-1) .* exp(-(1:floor(EEG.pnts/2)-1)/ed);
        as = [as(1) as 0 as(:,end:-1:1)];
        fc = as .* exp(1i*2*pi*rand(size(as)));
        noise = real(ifft(fc)) * EEG.pnts * 750; % scale noise
        
        EEG.data(chani,:,triali) = EEG.data(chani,:,triali) + noise;
    end
end

%% Plotting Section
if show_plots
    
    % Figure 1: Topographical maps of dipole projections
    figure(1), clf
    subplot(121)
    topoplotIndie(-lf.Gain(:,1,diploc1), EEG.chanlocs, 'numcontour', 0);
    title('Dipole 1 projection'), colorbar
    subplot(122)
    topoplotIndie(-lf.Gain(:,1,diploc2), EEG.chanlocs, 'numcontour', 0);
    title('Dipole 2 projection'), colorbar

    % --- FFT Analysis ---
    time_idx = (EEG.times >= time_range(1)) & (EEG.times <= time_range(2));
    n_points = sum(time_idx);
    frequencies = linspace(0, EEG.srate, n_points);

    % Figure 2: FFT spectrum for channels capturing each dipole
    [~, chan1] = max(abs(lf.Gain(:,1,diploc1)));
    [~, chan2] = max(abs(lf.Gain(:,1,diploc2)));
    
    fft_chan1 = fft(mean(EEG.data(chan1, time_idx, :), 3), n_points);
    fft_chan2 = fft(mean(EEG.data(chan2, time_idx, :), 3), n_points);
    
    freq_idx_plot = frequencies <= 50;
    pow_chan1 = 2*abs(fft_chan1(freq_idx_plot))/n_points;
    pow_chan2 = 2*abs(fft_chan2(freq_idx_plot))/n_points;
    
    figure(2), clf
    subplot(211)
    plot(frequencies(freq_idx_plot), pow_chan1, 'b', 'LineWidth', 2)
    hold on, plot([freq freq], [0 max(pow_chan1)*1.1], 'r--', 'LineWidth', 2)
    xlabel('Frequency (Hz)'), ylabel('Amplitude'), title(['FFT Spectrum - Channel ' num2str(chan1)]), xlim([0 50]), grid on
    
    subplot(212)
    plot(frequencies(freq_idx_plot), pow_chan2, 'r', 'LineWidth', 2)
    hold on, plot([freq freq], [0 max(pow_chan2)*1.1], 'r--', 'LineWidth', 2)
    xlabel('Frequency (Hz)'), ylabel('Amplitude'), title(['FFT Spectrum - Channel ' num2str(chan2)]), xlim([0 50]), grid on

    % Figure 3: Topoplot at specified frequency
    all_fft = abs(fft(mean(EEG.data(:, time_idx, :), 3), n_points, 2));
    [~, freq_idx] = min(abs(frequencies - freq));
    amp_at_freq = 2*all_fft(:, freq_idx)/n_points;
    
    [~, sorted_indices] = sort(amp_at_freq, 'descend');
    
    figure(3), clf
    topoplotIndie(amp_at_freq, EEG.chanlocs, 'numcontour', 0, 'electrodes', 'on');
    title(sprintf('Topography at %g Hz\nMax channels: %d and %d', freq, sorted_indices(1), sorted_indices(2)))
    colorbar
end

end