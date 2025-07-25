%% Inter-Site Phase Clustering (PLI) Tutorial
% This script demonstrates how to compute phase synchronization between EEG channels
% using the Inter-Site Phase Clustering method

%clear; close all;
addpath("G:\My Drive\Projects\signal_processing_mike_cohen\how_to\src")
% Define simulation parameters
sim_freq = 15;                       % Frequency in Hz
sim_phaselag = 0 * 2*pi;            % Phase lag in radians
dipole_1_loc = 101;                  % Index for dipole 1
dipole_2_loc = 201;                  % Index for dipole 2
activation_win = [0 2];              % Activation window in seconds
generate_plots = false;               % Generate plots (true) or not (false)

% Call the function to generate data
EEG = simulate_phaseLag_data(sim_freq, sim_phaselag, dipole_1_loc, dipole_2_loc, activation_win, generate_plots);


%% Create complex Morlet wavelet
target_freq = 15; % Hz (tagged frequency)
wavelet_time = -1:1/EEG.srate:1;
wavelet_cycles = 3;

% Gaussian width parameter
s = wavelet_cycles / (2*pi*target_freq);

% Complex Morlet wavelet
wavelet = exp(2i*pi*target_freq*wavelet_time) .* exp(-wavelet_time.^2/(2*s^2));

% Convolution parameters
nwave = length(wavelet);
ndata = EEG.pnts;
nconv = nwave + ndata - 1;
half_wave = floor(nwave/2);

% Normalize wavelet
waveletX = fft(wavelet, nconv);
waveletX = waveletX / max(waveletX);

%% Extract phase data for all channels
phase_data = zeros(EEG.nbchan, EEG.pnts, EEG.trials);

for chani = 1:EEG.nbchan
    for triali = 1:EEG.trials
        % FFT of data
        dataX = fft(EEG.data(chani, :, triali), nconv);
        
        % Convolution
        conv_result = ifft(waveletX .* dataX, nconv);
        conv_result = conv_result(half_wave+1:end-half_wave);
        
        % Extract phase
        phase_data(chani, :, triali) = angle(conv_result);
    end
end

%% Compute PLI between all channel pairs

% define time windows
tidx = dsearchn(EEG.times',[0 2]');

PLI_matrix = zeros(EEG.nbchan, EEG.nbchan);

for chan1 = 1:EEG.nbchan
    for chan2 = 1:EEG.nbchan
        % Phase differences across all trials and time points
        phase_diff = phase_data(chan1, tidx(1):tidx(2), :) - phase_data(chan2, tidx(1):tidx(2), :);
        phase_diff = phase_diff(:);

        % PLI calculation
        complex_diff = exp(1i * phase_diff);
        PLI_matrix(chan1, chan2) = abs(mean(sign(imag(complex_diff))));
    end
end

% Set diagonal to 1 (perfect synchronization with itself)
%PLI_matrix = PLI_matrix + eye(EEG.nbchan);

%% Visualize PLI matrix
figure('Position', [100, 100, 600, 500]);

imagesc(PLI_matrix);
colorbar;
%colormap('hot');
caxis([0 0.05]);
xlabel('Channel');
ylabel('Channel');
title('Inter-Site Phase Clustering (PLI) Matrix');
axis square;

% Add grid for better visualization
hold on;
for i = 0.5:1:EEG.nbchan+0.5
    plot([i i], [0.5 EEG.nbchan+0.5], 'k-', 'LineWidth', 0.5);
    plot([0.5 EEG.nbchan+0.5], [i i], 'k-', 'LineWidth', 0.5);
end

