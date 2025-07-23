%% ========================================================================
%% TUTORIAL: Principal Component Analysis (PCA) for EEG Data
%% ========================================================================
%
% This tutorial demonstrates how to apply PCA to EEG data to understand
% phase-locked vs. total variance in neural signals.
%
% WHAT IS PCA?
% Principal Component Analysis (PCA) is a dimensionality reduction technique
% that finds the directions of maximum variance in high-dimensional data.
% 
% KEY CONCEPTS:
% 1. PCA transforms correlated variables into uncorrelated components
% 2. Components are ordered by variance explained (eigenvalues)
% 3. First component captures most variance, second captures second most, etc.
% 4. Covariance matrix contains information about how channels vary together
% 5. Eigenvectors of covariance matrix are the principal components
% 6. Eigenvalues tell us how much variance each component explains
%
% WHY USE PCA WITH EEG?
% - EEG has many channels that are often correlated
% - PCA helps identify dominant patterns of brain activity
% - Reduces data complexity while preserving important information
%
% TWO TYPES OF PCA IN EEG:
% 1. PHASE-LOCKED PCA: Applied to trial-averaged data (ERP)
%    - Shows consistent activity across trials
%    - Captures evoked responses
% 
% 2. TOTAL PCA: Applied to single-trial covariances
%    - Shows all variance (evoked + induced + noise)
%    - Captures both consistent and variable activity
%
%% ========================================================================

% Clear workspace for clean start
close all; clear; clc

fprintf('=== PCA TUTORIAL FOR EEG DATA ===\n\n')

%% STEP 1: Load and Prepare EEG Data
fprintf('STEP 1: Loading EEG data...\n')

load sampleEEGdata.mat
EEG.data = double(EEG.data); % Ensure data is in double precision

% Display data structure
fprintf('Data dimensions: %d channels x %d time points x %d trials\n', ...
    EEG.nbchan, EEG.pnts, EEG.trials)
fprintf('Time range: %.0f to %.0f ms\n', EEG.times(1), EEG.times(end))
fprintf('Sampling rate: %.0f Hz\n\n', EEG.srate)

%% STEP 2: Define Analysis Time Window
fprintf('STEP 2: Defining analysis time window...\n')

% Choose time window for analysis (0 to 800 ms post-stimulus)
time_window = [0 800]; % in milliseconds
tidx = dsearchn(EEG.times', time_window');

fprintf('Analyzing time window: %.0f to %.0f ms\n', ...
    EEG.times(tidx(1)), EEG.times(tidx(2)))
fprintf('Time indices: %d to %d\n\n', tidx(1), tidx(2))

%% STEP 3: Compute Phase-Locked (ERP) Covariance Matrix
fprintf('STEP 3: Computing phase-locked covariance matrix...\n')

% Phase-locked analysis: average across trials first, then compute covariance
% This captures only activity that is consistent across trials (ERPs)

% Step 3a: Compute trial average (ERP)
erp = mean(EEG.data(:, tidx(1):tidx(2), :), 3);
fprintf('ERP computed: averaged across %d trials\n', EEG.trials)

% Step 3b: Remove channel means (centering the data)
data_centered = erp - mean(erp, 2);
fprintf('Data centered: removed channel means\n')

% Step 3c: Compute covariance matrix
% Covariance shows how channels co-vary together
cov_phaselocked = (data_centered * data_centered') / (size(data_centered,2) - 1);
fprintf('Phase-locked covariance matrix computed (%dx%d)\n\n', ...
    size(cov_phaselocked,1), size(cov_phaselocked,2))

%% STEP 4: Compute Total (Single-Trial) Covariance Matrix
fprintf('STEP 4: Computing total covariance matrix...\n')

% Total analysis: compute covariance for each trial, then average
% This captures all variance (evoked + induced + noise)

cov_total = zeros(EEG.nbchan, EEG.nbchan);

fprintf('Computing covariance for each trial:\n')
for triali = 1:EEG.trials
    % Extract data from this trial and compute covariance
    trial_data = EEG.data(:, tidx(1):tidx(2), triali)';
    trial_cov = cov(trial_data);
    
    % Add to running sum
    cov_total = cov_total + trial_cov;
    
    % Progress indicator
    if mod(triali, 20) == 0
        fprintf('  Processed %d/%d trials\n', triali, EEG.trials)
    end
end

% Average across trials
cov_total = cov_total / EEG.trials;
fprintf('Total covariance matrix computed (averaged across %d trials)\n\n', EEG.trials)

%% STEP 5: Visualize Covariance Matrices
fprintf('STEP 5: Visualizing covariance matrices...\n')

figure(1), clf
set(gcf, 'Position', [100 100 800 350])

% Phase-locked covariance
subplot(1,2,1)
imagesc(cov_phaselocked)
axis square
colorbar
set(gca, 'clim', [-1 1]*2)
xlabel('Channels')
ylabel('Channels')
title('Phase-Locked Covariance (ERP)')
colormap(jet)

% Total covariance
subplot(1,2,2)
imagesc(cov_total)
axis square
colorbar
set(gca, 'clim', [-1 1]*100)
xlabel('Channels')
ylabel('Channels')
title('Total Covariance (Single Trials)')

fprintf('Covariance matrices displayed\n')
fprintf('Notice: Total covariance has much larger values (more variance)\n\n')

%% STEP 6: Perform Phase-Locked PCA
fprintf('STEP 6: Performing phase-locked PCA...\n')

% PCA via eigendecomposition of covariance matrix
[V_PL, L_PL] = eig(cov_phaselocked);

% Sort eigenvalues and eigenvectors in descending order
[L_PL, sidx] = sort(diag(L_PL), 'descend');
V_PL = V_PL(:, sidx);

% Convert eigenvalues to percentage of total variance
L_PL_pct = 100 * L_PL / sum(L_PL);

fprintf('Phase-locked PCA completed:\n')
fprintf('  First component explains %.1f%% of variance\n', L_PL_pct(1))
fprintf('  First 5 components explain %.1f%% of variance\n', sum(L_PL_pct(1:5)))

% Compute principal component time series
% Apply the first component filter to the trial-averaged data
compts_PL = V_PL(:,1)' * mean(EEG.data, 3);

%% STEP 7: Perform Total PCA
fprintf('\nSTEP 7: Performing total PCA...\n')

% PCA of total covariance matrix
[V_TT, L_TT] = eig(cov_total);
[L_TT, sidx] = sort(diag(L_TT), 'descend');
V_TT = V_TT(:, sidx);
L_TT_pct = 100 * L_TT / sum(L_TT);

fprintf('Total PCA completed:\n')
fprintf('  First component explains %.1f%% of variance\n', L_TT_pct(1))
fprintf('  First 5 components explain %.1f%% of variance\n', sum(L_TT_pct(1:5)))

% Compute principal component time series for total PCA
% Need to reshape 3D data to 2D, apply filter, then reshape back
data2d = reshape(EEG.data, EEG.nbchan, []);
comptmp = V_TT(:,1)' * data2d;
compts_TT = reshape(comptmp, [EEG.pnts, EEG.trials]);

% Compute trial average of the component
compTT_erp = mean(compts_TT, 2);

%% STEP 8: Visualize PCA Results
fprintf('\nSTEP 8: Visualizing PCA results...\n')

figure(2), clf
set(gcf, 'Position', [200 50 900 700])

% Scree plot - variance explained by each component
subplot(3,2,1)
plot(1:20, L_PL_pct(1:20), 'ks-', 'markerfacecolor', 'b', 'markersize', 8)
hold on
plot(1:20, L_TT_pct(1:20), 'ro-', 'markerfacecolor', 'r', 'markersize', 8)
xlabel('Component Number')
ylabel('Variance Explained (%)')
title('Scree Plot: Variance per Component')
legend({'Phase-locked', 'Total'}, 'Location', 'northeast')
set(gca, 'xlim', [0.5 20])
grid on

% Cumulative variance plot
subplot(3,2,2)
plot(1:20, cumsum(L_PL_pct(1:20)), 'ks-', 'markerfacecolor', 'b', 'markersize', 8)
hold on
plot(1:20, cumsum(L_TT_pct(1:20)), 'ro-', 'markerfacecolor', 'r', 'markersize', 8)
xlabel('Component Number')
ylabel('Cumulative Variance (%)')
title('Cumulative Variance Explained')
legend({'Phase-locked', 'Total'}, 'Location', 'southeast')
set(gca, 'xlim', [0.5 20])
grid on

% Topographical maps of first components
subplot(3,2,3)
topoplotIndie(V_PL(:,1), EEG.chanlocs, 'numcontour', 0, 'shading', 'interp', 'electrodes', 'off');
title('Phase-Locked PC1 Topography')
colorbar

subplot(3,2,4)
topoplotIndie(V_TT(:,1), EEG.chanlocs, 'numcontour', 0, 'shading', 'interp', 'electrodes', 'off');
title('Total PC1 Topography')
colorbar

% Time courses of first components
subplot(3,1,3)
plot(EEG.times, compts_PL, 'b-', 'linewidth', 2)
hold on
plot(EEG.times, compTT_erp, 'r-', 'linewidth', 2)
xlabel('Time (ms)')
ylabel('Component Activity (μV)')
title('First Principal Component Time Courses')
legend({'Phase-locked PC1', 'Total PC1'}, 'Location', 'best')
set(gca, 'xlim', [-200 1300])
grid on

%.