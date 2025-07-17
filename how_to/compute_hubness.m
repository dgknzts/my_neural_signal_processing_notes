%% Connectivity Hubs and Graph Theory Analysis
% This script demonstrates connectivity analysis using ISPC and graph theory metrics
% Based on neural signal processing principles

clear; close all;
addpath("G:\My Drive\Projects\signal_processing_mike_cohen\how_to\src")

%% Generate simulated EEG data
sim_freq = 15;
sim_phaselag = 0.25 * 2*pi;
dipole_1_loc = 32;  % PO7
dipole_2_loc = 86;  % PO8
activation_win = [0 2];
generate_plots = false;

EEG = simulate_phaseLag_data(sim_freq, sim_phaselag, dipole_1_loc, dipole_2_loc, ...
    activation_win, 'noise_level', 100, 'signal_amp', [2, 2], 'show_plots', generate_plots);

%% Connectivity Analysis Parameters
target_freq = 15;
wavelet_time = -0.7:1/EEG.srate:0.7;  % 1.4 seconds total
wavelet_cycles = 12;
time_window = [0 2];

% Gaussian width parameter
s = wavelet_cycles / (2*pi*target_freq);

% Complex Morlet wavelet
wavelet = exp(2i*pi*target_freq*wavelet_time) .* exp(-wavelet_time.^2/(2*s^2));

figure;
subplot(211); plot(wavelet_time, real(wavelet)); title('Time domain');
subplot(212); plot(abs(fft(wavelet))); title('Frequency domain');

% Convolution parameters
nwave = length(wavelet);
ndata = EEG.pnts;
nconv = nwave + ndata - 1;
half_wave = floor(nwave/2);

% Normalize wavelet
waveletX = fft(wavelet, nconv);
waveletX = waveletX / max(waveletX);

%% Extract phase data for all channels
fprintf('Extracting phase data...\n');
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

%% Compute connectivity matrices
fprintf('Computing connectivity matrices...\n');
tidx = dsearchn(EEG.times', time_window');

% Initialize matrices
ispc_matrix = zeros(EEG.nbchan, EEG.nbchan);
pli_matrix = zeros(EEG.nbchan, EEG.nbchan);

for chan1 = 1:EEG.nbchan
    for chan2 = chan1+1:EEG.nbchan
        % Phase differences across trials and time
        phase_diff = phase_data(chan1, tidx(1):tidx(2), :) - phase_data(chan2, tidx(1):tidx(2), :);
        phase_diff = phase_diff(:);
        
        % ISPC calculation
        complex_diff = exp(1i * phase_diff);
        ispc_val = abs(mean(complex_diff));
        
        % PLI calculation
        pli_val = abs(mean(sign(imag(complex_diff))));
        
        % Fill symmetric matrices
        ispc_matrix(chan1, chan2) = ispc_val;
        ispc_matrix(chan2, chan1) = ispc_val;
        pli_matrix(chan1, chan2) = pli_val;
        pli_matrix(chan2, chan1) = pli_val;
    end
end

%% Graph Theory Analysis
fprintf('Computing graph theory metrics...\n');

% Choose connectivity matrix for analysis
conn_matrix = ispc_matrix;  % or pli_matrix

% Thresholding methods
unique_vals = nonzeros(triu(conn_matrix));
thresh_mean = mean(unique_vals);
thresh_std = thresh_mean + std(unique_vals);
thresh_percentile = prctile(unique_vals, 95);

% Apply threshold
thresh_method = 'std';  % 'mean', 'std', or 'percentile'
switch thresh_method
    case 'mean'
        threshold = thresh_mean;
    case 'std'
        threshold = thresh_std;
    case 'percentile'
        threshold = thresh_percentile;
end
%threshold = 0.078;
% Create binary adjacency matrix
adj_matrix = conn_matrix > threshold;

% Graph theory metrics
hubness = sum(adj_matrix, 2) / (EEG.nbchan - 1);  % Degree centrality
clustering = compute_clustering(adj_matrix);
path_length = compute_path_length(adj_matrix);

% Global metrics
global_clustering = mean(clustering);
global_path_length = mean(path_length(~isinf(path_length)));
density = sum(adj_matrix(:)) / (EEG.nbchan * (EEG.nbchan - 1));

%% Visualizations
fprintf('Creating visualizations...\n');

% Figure 1: Connectivity matrices
figure('Position', [100, 100, 1200, 400]);

subplot(1,3,1)
imagesc(ispc_matrix);
colorbar;
title('ISPC Matrix');
xlabel('Channel'); ylabel('Channel');
axis square;

subplot(1,3,2)
imagesc(pli_matrix);
colorbar;
title('PLI Matrix');
xlabel('Channel'); ylabel('Channel');
axis square;

subplot(1,3,3)
imagesc(adj_matrix);
colorbar;
title(sprintf('Thresholded Matrix (%.3f)', threshold));
xlabel('Channel'); ylabel('Channel');
axis square;

% Figure 2: Threshold analysis
figure('Position', [100, 200, 1200, 400]);

subplot(1,3,1)
histogram(unique_vals, 50);
hold on;
plot([threshold threshold], ylim, 'r--', 'LineWidth', 2);
xlabel('Connectivity Strength');
ylabel('Count');
title('Connectivity Distribution');
legend('Distribution', 'Threshold');

subplot(1,3,2)
topoplotIndie(hubness, EEG.chanlocs, 'numcontour', 0);
title('Degree Centrality (Hubness)');
colorbar;

subplot(1,3,3)
topoplotIndie(clustering, EEG.chanlocs, 'numcontour', 0);
title('Clustering Coefficient');
colorbar;

% Figure 3: Graph metrics comparison
figure('Position', [100, 300, 800, 600]);

subplot(2,2,1)
bar(hubness);
xlabel('Channel');
ylabel('Degree Centrality');
title('Hub Strength per Channel');

subplot(2,2,2)
bar(clustering);
xlabel('Channel');
ylabel('Clustering Coefficient');
title('Local Clustering per Channel');

subplot(2,2,3)
scatter(hubness, clustering);
xlabel('Degree Centrality');
ylabel('Clustering Coefficient');
title('Hub vs Clustering Relationship');

subplot(2,2,4)
% Network statistics
stats_text = {
    sprintf('Density: %.3f', density),
    sprintf('Global Clustering: %.3f', global_clustering),
    sprintf('Global Path Length: %.3f', global_path_length),
    sprintf('Threshold: %.3f', threshold),
    sprintf('Connected Edges: %d', sum(adj_matrix(:))/2)
};
text(0.1, 0.8, stats_text, 'FontSize', 12, 'Units', 'normalized');
title('Network Statistics');
axis off;

% Figure 4: Top hubs analysis
[~, hub_idx] = sort(hubness, 'descend');
top_hubs = hub_idx(1:5);

figure('Position', [100, 400, 1000, 300]);
for i = 1:5
    subplot(1,5,i);
    hub_pattern = zeros(EEG.nbchan, 1);
    hub_pattern(adj_matrix(top_hubs(i), :)) = 1;
    topoplotIndie(hub_pattern, EEG.chanlocs, 'numcontour', 0);
    title(sprintf('Hub %d (Ch %d)', i, top_hubs(i)));
end
sgtitle('Top 5 Connectivity Hubs');

%% Results Summary
fprintf('\n=== CONNECTIVITY ANALYSIS RESULTS ===\n');
fprintf('Frequency analyzed: %.1f Hz\n', target_freq);
fprintf('Time window: [%.1f, %.1f] s\n', time_window);
fprintf('Threshold method: %s (%.3f)\n', thresh_method, threshold);
fprintf('Network density: %.3f\n', density);
fprintf('Global clustering: %.3f\n', global_clustering);
fprintf('Global path length: %.3f\n', global_path_length);
fprintf('Number of connected edges: %d\n', sum(adj_matrix(:))/2);

fprintf('\nTop 5 hubs (by degree centrality):\n');
for i = 1:5
    fprintf('  %d. Channel %d: %.3f\n', i, top_hubs(i), hubness(top_hubs(i)));
end

%% Helper Functions
function clustering = compute_clustering(adj_matrix)
    % Compute clustering coefficient for each node
    n = size(adj_matrix, 1);
    clustering = zeros(n, 1);
    
    for i = 1:n
        neighbors = find(adj_matrix(i, :));
        ki = length(neighbors);
        
        if ki < 2
            clustering(i) = 0;
        else
            % Count triangles
            subgraph = adj_matrix(neighbors, neighbors);
            triangles = sum(subgraph(:)) / 2;
            clustering(i) = 2 * triangles / (ki * (ki - 1));
        end
    end
end

function path_length = compute_path_length(adj_matrix)
    % Compute average shortest path length for each node
    n = size(adj_matrix, 1);
    path_length = zeros(n, 1);
    
    % Convert to distance matrix (1 for connected, Inf for disconnected)
    dist_matrix = double(adj_matrix);
    dist_matrix(dist_matrix == 0) = Inf;
    dist_matrix(1:n+1:end) = 0;  % diagonal = 0
    
    % Floyd-Warshall algorithm
    for k = 1:n
        for i = 1:n
            for j = 1:n
                if dist_matrix(i,k) + dist_matrix(k,j) < dist_matrix(i,j)
                    dist_matrix(i,j) = dist_matrix(i,k) + dist_matrix(k,j);
                end
            end
        end
    end
    
    % Compute average path length for each node
    for i = 1:n
        paths = dist_matrix(i, :);
        paths = paths(paths ~= 0 & ~isinf(paths));  % exclude self and disconnected
        if ~isempty(paths)
            path_length(i) = mean(paths);
        else
            path_length(i) = Inf;
        end
    end
end