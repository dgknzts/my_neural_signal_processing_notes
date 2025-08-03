%% Simple EEG Spatial Filtering: ZCA + PCA + GED
% Basic tutorial showing the essential steps

close all; clear; clc

%% Load data
load sampleEEGdata.mat
EEG.data = double(EEG.data);

fprintf('Data: %d channels, %d timepoints, %d trials\n', EEG.nbchan, EEG.pnts, EEG.trials);

%% Save original data for comparison
original_data = EEG.data;

%% STEP 1: ZCA Whitening (remove channel correlations)

% Get baseline period for noise estimation
baseline_window = [-1000 0]; % ms
tidx = dsearchn(EEG.times', baseline_window');

% Extract baseline data from all trials
baseline_data = reshape(EEG.data(:, tidx(1):tidx(2), :), EEG.nbchan, []);
% Remove mean from each channel
baseline_data = bsxfun(@minus, baseline_data, mean(baseline_data, 2));

% Compute noise covariance matrix
noise_cov = baseline_data * baseline_data' / size(baseline_data, 2);

% Eigendecomposition of noise
[V_noise, D_noise] = eig(noise_cov);
eigenvals = diag(D_noise);
eigenvals_pct = 100 * eigenvals / sum(eigenvals);

% Keep components with >0.1% variance
keep_idx = eigenvals_pct > 0.1;
fprintf('ZCA: keeping %d/%d components\n', sum(keep_idx), EEG.nbchan);

% Create whitening matrix
whitening_matrix = V_noise(:, keep_idx) * diag(eigenvals(keep_idx).^(-0.5)) * V_noise(:, keep_idx)';

% Apply ZCA whitening
EEG.data = reshape(whitening_matrix * reshape(EEG.data, EEG.nbchan, []), size(EEG.data));

%% STEP 2: PCA Compression (reduce dimensions)

% Compute covariance across all whitened data
data_cov = zeros(EEG.nbchan);
for trial = 1:EEG.trials
    % Remove linear trend from each trial
    trial_data = detrend(EEG.data(:, :, trial)', 'linear')';
    data_cov = data_cov + trial_data * trial_data' / EEG.pnts;
end
data_cov = data_cov / EEG.trials;

% PCA decomposition
[V_pca, D_pca] = eig(data_cov);
[eigenvals_pca, sort_idx] = sort(diag(D_pca), 'descend');
V_pca = V_pca(:, sort_idx);

% Keep components with >1% variance
eigenvals_pct = 100 * eigenvals_pca / sum(eigenvals_pca);
pca_keep = eigenvals_pct > 1;
n_comps = sum(pca_keep);
fprintf('PCA: keeping %d components (%.1f%% variance)\n', n_comps, sum(eigenvals_pct(1:n_comps)));

% Project data to PCA space
pca_weights = V_pca(:, pca_keep)';
compressed_data = reshape(pca_weights * reshape(EEG.data, EEG.nbchan, []), [n_comps, EEG.pnts, EEG.trials]);

%% STEP 3: GED Spatial Filtering (maximize SNR)

% Create filtered data at target frequency
target_freq = 11; % Hz
filter_width = 1; % Hz
filtered_data = filterFGx(compressed_data, EEG.srate, target_freq, filter_width);

% Compute covariance matrices
broadband_cov = zeros(n_comps);
signal_cov = zeros(n_comps);

for trial = 1:EEG.trials
    % Broadband covariance (noise)
    trial_broad = detrend(compressed_data(:, :, trial)', 'linear')';
    broadband_cov = broadband_cov + trial_broad * trial_broad' / EEG.pnts;
    
    % Filtered covariance (signal)
    trial_filt = detrend(filtered_data(:, :, trial)', 'linear')';
    signal_cov = signal_cov + trial_filt * trial_filt' / EEG.pnts;
end

broadband_cov = broadband_cov / EEG.trials;
signal_cov = signal_cov / EEG.trials;

% Solve generalized eigenvalue problem
[spatial_filters, lambda] = eig(signal_cov, broadband_cov);
[lambda, sort_idx] = sort(diag(lambda), 'descend');
spatial_filters = spatial_filters(:, sort_idx);

fprintf('Best component SNR: %.2f\n', lambda(1));

% Extract time series of best component
best_filter = spatial_filters(:, 1);
component_ts = reshape(reshape(compressed_data, n_comps, [])' * best_filter, EEG.pnts, EEG.trials);

% Compute topography (forward model)
topography = best_filter' * signal_cov * pca_weights;

% Fix sign for interpretability
[~, max_idx] = max(abs(topography));
if topography(max_idx) < 0
    topography = -topography;
    component_ts = -component_ts;
end

%% Create comparison: Original vs ZCA preprocessed

% Run same analysis on original data (without ZCA)
fprintf('\nRunning analysis without ZCA for comparison...\n');

% PCA on original data
orig_cov = zeros(EEG.nbchan);
for trial = 1:EEG.trials
    trial_data = detrend(original_data(:, :, trial)', 'linear')';
    orig_cov = orig_cov + trial_data * trial_data' / EEG.pnts;
end
orig_cov = orig_cov / EEG.trials;

[V_orig, D_orig] = eig(orig_cov);
[evals_orig, sort_idx] = sort(diag(D_orig), 'descend');
V_orig = V_orig(:, sort_idx);

% Keep same number of components for fair comparison
orig_weights = V_orig(:, 1:n_comps)';
orig_compressed = reshape(orig_weights * reshape(original_data, EEG.nbchan, []), [n_comps, EEG.pnts, EEG.trials]);

% Filter original compressed data
orig_filtered = filterFGx(orig_compressed, EEG.srate, target_freq, filter_width);

% GED on original data
orig_broad_cov = zeros(n_comps);
orig_signal_cov = zeros(n_comps);

for trial = 1:EEG.trials
    trial_broad = detrend(orig_compressed(:, :, trial)', 'linear')';
    orig_broad_cov = orig_broad_cov + trial_broad * trial_broad' / EEG.pnts;
    
    trial_filt = detrend(orig_filtered(:, :, trial)', 'linear')';
    orig_signal_cov = orig_signal_cov + trial_filt * trial_filt' / EEG.pnts;
end

orig_broad_cov = orig_broad_cov / EEG.trials;
orig_signal_cov = orig_signal_cov / EEG.trials;

[orig_filters, orig_lambda] = eig(orig_signal_cov, orig_broad_cov);
[orig_lambda, sort_idx] = sort(diag(orig_lambda), 'descend');
orig_filters = orig_filters(:, sort_idx);

% Extract original component
orig_best_filter = orig_filters(:, 1);
orig_component_ts = reshape(reshape(orig_compressed, n_comps, [])' * orig_best_filter, EEG.pnts, EEG.trials);
orig_topography = orig_best_filter' * orig_signal_cov * orig_weights;

% Fix sign
[~, max_idx] = max(abs(orig_topography));
if orig_topography(max_idx) < 0
    orig_topography = -orig_topography;
    orig_component_ts = -orig_component_ts;
end

fprintf('Original data SNR: %.2f\n', orig_lambda(1));

%% Plot results
figure('Position', [100, 100, 1200, 400]);

% Topography without ZCA
subplot(1, 3, 1);
topoplotIndie(orig_topography, EEG.chanlocs, 'numcontour', 0);
title(sprintf('Without ZCA\nSNR = %.2f', orig_lambda(1)));
colorbar;

% Topography with ZCA
subplot(1, 3, 2);
topoplotIndie(topography, EEG.chanlocs, 'numcontour', 0);
title(sprintf('With ZCA\nSNR = %.2f', lambda(1)));
colorbar;

% Component ERPs comparison
subplot(1, 3, 3);
plot(EEG.times, mean(orig_component_ts, 2), 'b-', 'LineWidth', 2, 'DisplayName', 'Without ZCA');
hold on;
plot(EEG.times, mean(component_ts, 2), 'r-', 'LineWidth', 2, 'DisplayName', 'With ZCA');
xlabel('Time (ms)');
ylabel('Component Amplitude');
title('Component ERPs');
xlim([-200, 1000]);
legend('Location', 'best');
grid on;

sgtitle(sprintf('ZCA + PCA + GED Results (Target: %.0f Hz)', target_freq), 'FontSize', 14);
