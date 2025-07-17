function EEG = simulate_phaseLag_data(freq, phaselag_rad, diploc1, diploc2, time_range, varargin)
% Enhanced simulation of EEG data from two phase-lagged dipoles
%
% Inputs:
%   freq          - Signal frequency in Hz
%   phaselag_rad  - Phase difference between dipoles in radians
%   diploc1       - Index of the first dipole location
%   diploc2       - Index of the second dipole location  
%   time_range    - 1x2 vector for signal activation window [start, end] in seconds
%
% Optional parameters (name-value pairs):
%   'show_plots'     - Boolean to control plotting (default: false)
%   'noise_level'    - Background noise scaling factor (default: 750)
%   'pink_exp'       - Pink noise exponential decay (default: 50)
%   'signal_amp'     - Signal amplitude [amp1, amp2] (default: [6, 6])
%   'amp_jitter'     - Trial-to-trial amplitude variability (default: 0.05)
%   'n_trials'       - Number of trials (default: 80)
%   'window_type'    - Activation window type: 'hann', 'rect', 'gauss' (default: 'hann')
%   'smooth_duration'- Smoothing window duration in seconds (default: 0.1)
%   'freq_jitter'    - Frequency jitter std (default: 0)
%   'phase_jitter'   - Phase jitter std in radians (default: 0)

% Parse optional inputs
p = inputParser;
addParameter(p, 'show_plots', false, @islogical);
addParameter(p, 'noise_level', 750, @isnumeric);
addParameter(p, 'pink_exp', 50, @isnumeric);
addParameter(p, 'signal_amp', [6, 6], @isnumeric);
addParameter(p, 'amp_jitter', 0.05, @isnumeric);
addParameter(p, 'n_trials', 80, @isnumeric);
addParameter(p, 'window_type', 'hann', @ischar);
addParameter(p, 'smooth_duration', 0.1, @isnumeric);
addParameter(p, 'freq_jitter', 0, @isnumeric);
addParameter(p, 'phase_jitter', 0, @isnumeric);
parse(p, varargin{:});

% Extract parameters
show_plots = p.Results.show_plots;
noise_level = p.Results.noise_level;
pink_exp = p.Results.pink_exp;
signal_amp = p.Results.signal_amp;
amp_jitter = p.Results.amp_jitter;
n_trials = p.Results.n_trials;
window_type = p.Results.window_type;
smooth_duration = p.Results.smooth_duration;
freq_jitter = p.Results.freq_jitter;
phase_jitter = p.Results.phase_jitter;

% Ensure signal_amp has two values
if length(signal_amp) == 1
    signal_amp = [signal_amp, signal_amp];
end

%% Setup
addpath("src");
load emptyEEG.mat

% EEG parameters
EEG.times = (-1:1/EEG.srate:2.4);
EEG.pnts = numel(EEG.times);
EEG.trials = n_trials;
EEG.data = zeros(EEG.nbchan, EEG.pnts, EEG.trials);

% Create activation window
rect_win = (EEG.times >= time_range(1)) & (EEG.times <= time_range(2));

% Apply different window types
switch lower(window_type)
    case 'rect'
        smooth_len = round(smooth_duration * EEG.srate);
        activation_win = conv(double(rect_win), hann(smooth_len)/sum(hann(smooth_len)), 'same');
    case 'hann'
        win_samples = sum(rect_win);
        full_win = zeros(size(EEG.times));
        full_win(rect_win) = hann(win_samples);
        activation_win = full_win;
    case 'gauss'
        win_center = mean(time_range);
        win_width = diff(time_range) / 4; % 4-sigma window
        activation_win = exp(-((EEG.times - win_center).^2) / (2 * win_width^2));
    otherwise
        error('Unknown window type: %s', window_type);
end

%% Data Simulation
rng(1); % Reproducibility
for triali = 1:EEG.trials
    
    % Background dipole activity (spatially uncorrelated noise)
    dipole_data = 0.02 * randn(size(lf.Gain, 3), EEG.pnts);
    
    % Trial-specific parameters with jitter
    trial_freq1 = freq + freq_jitter * randn;
    trial_freq2 = freq + freq_jitter * randn;
    trial_phase1 = phase_jitter * randn;
    trial_phase2 = phaselag_rad + phase_jitter * randn;
    trial_amp1 = signal_amp(1) + amp_jitter * randn;
    trial_amp2 = signal_amp(2) + amp_jitter * randn;
    
    % Generate signals
    sig1 = trial_amp1 * sin(2*pi*trial_freq1*EEG.times + trial_phase1) .* activation_win;
    sig2 = trial_amp2 * sin(2*pi*trial_freq2*EEG.times + trial_phase2) .* activation_win;
    
    % Assign to dipoles
    dipole_data(diploc1, :) = dipole_data(diploc1, :) + sig1;
    dipole_data(diploc2, :) = dipole_data(diploc2, :) + sig2;
    
    % Project to scalp
    EEG.data(:, :, triali) = squeeze(lf.Gain(:, 1, :)) * dipole_data;
end

%% Add realistic noise
if noise_level > 0
    EEG.data = add_realistic_noise(EEG.data, EEG.srate, noise_level, pink_exp);
end

%% Plotting
if show_plots
    create_diagnostic_plots(EEG, lf, diploc1, diploc2, freq, time_range);
end

end

%% Helper Functions
function data_with_noise = add_realistic_noise(data, srate, noise_level, pink_exp)
% Add realistic pink noise to EEG data
[n_chan, n_pnts, n_trials] = size(data);
data_with_noise = data;

for chani = 1:n_chan
    for triali = 1:n_trials
        % Generate pink noise in frequency domain
        half_pnts = floor(n_pnts/2);
        freqs = 1:half_pnts-1;
        
        % Power law decay
        as = rand(1, length(freqs)) .* exp(-freqs/pink_exp);
        as = [as(1), as, 0, as(end:-1:1)];
        
        % Add random phases
        fc = as .* exp(1i*2*pi*rand(size(as)));
        
        % Convert to time domain
        noise = real(ifft(fc)) * n_pnts * noise_level;
        
        data_with_noise(chani, :, triali) = data_with_noise(chani, :, triali) + noise;
    end
end
end

function create_diagnostic_plots(EEG, lf, diploc1, diploc2, freq, time_range)
% Create comprehensive diagnostic plots

% Figure 1: Dipole topographies
figure(1), clf
subplot(121)
topoplotIndie(-lf.Gain(:,1,diploc1), EEG.chanlocs, 'numcontour', 0);
title('Dipole 1 projection'), colorbar
subplot(122)
topoplotIndie(-lf.Gain(:,1,diploc2), EEG.chanlocs, 'numcontour', 0);
title('Dipole 2 projection'), colorbar

% Figure 2: Time series from peak channels
[~, chan1] = max(abs(lf.Gain(:,1,diploc1)));
[~, chan2] = max(abs(lf.Gain(:,1,diploc2)));

figure(2), clf
subplot(211)
plot(EEG.times, squeeze(mean(EEG.data(chan1, :, :), 3)), 'LineWidth', 2)
title(sprintf('Channel %d (Dipole 1 peak)', chan1))
xlabel('Time (s)'), ylabel('Amplitude')
grid on

subplot(212)
plot(EEG.times, squeeze(mean(EEG.data(chan2, :, :), 3)), 'LineWidth', 2)
title(sprintf('Channel %d (Dipole 2 peak)', chan2))
xlabel('Time (s)'), ylabel('Amplitude')
grid on

% Figure 3: Frequency analysis
time_idx = (EEG.times >= time_range(1)) & (EEG.times <= time_range(2));
n_points = sum(time_idx);
frequencies = linspace(0, EEG.srate, n_points);

fft_chan1 = fft(mean(EEG.data(chan1, time_idx, :), 3), n_points);
fft_chan2 = fft(mean(EEG.data(chan2, time_idx, :), 3), n_points);

freq_idx_plot = frequencies <= 50;
pow_chan1 = 2*abs(fft_chan1(freq_idx_plot))/n_points;
pow_chan2 = 2*abs(fft_chan2(freq_idx_plot))/n_points;

figure(3), clf
subplot(211)
plot(frequencies(freq_idx_plot), pow_chan1, 'b', 'LineWidth', 2)
hold on
plot([freq freq], [0 max(pow_chan1)*1.1], 'r--', 'LineWidth', 2)
xlabel('Frequency (Hz)'), ylabel('Amplitude')
title(sprintf('FFT - Channel %d', chan1))
xlim([0 50]), grid on

subplot(212)
plot(frequencies(freq_idx_plot), pow_chan2, 'r', 'LineWidth', 2)
hold on
plot([freq freq], [0 max(pow_chan2)*1.1], 'r--', 'LineWidth', 2)
xlabel('Frequency (Hz)'), ylabel('Amplitude')
title(sprintf('FFT - Channel %d', chan2))
xlim([0 50]), grid on

% Figure 4: Topography at target frequency
all_fft = abs(fft(mean(EEG.data(:, time_idx, :), 3), n_points, 2));
[~, freq_idx] = min(abs(frequencies - freq));
amp_at_freq = 2*all_fft(:, freq_idx)/n_points;

figure(4), clf
topoplotIndie(amp_at_freq, EEG.chanlocs, 'numcontour', 0, 'electrodes', 'on');
title(sprintf('Topography at %g Hz', freq))
colorbar
end