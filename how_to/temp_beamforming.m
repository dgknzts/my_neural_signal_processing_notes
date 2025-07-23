%% Source Reconstruction and Synchronization Analysis
% This script performs source reconstruction using beamforming (LCMV) 
% and analyzes synchronization between reconstructed dipole sources

clear; close all;
addpath("src")

load emptyEEG % Load leadfield
%% Generate simulated EEG data
sim_freq = 15;
sim_phaselag = 0.25 * 2*pi;  % Try non-zero phase lag first
dipole_1_loc = 56;
dipole_2_loc = 32;
activation_win = [0 2];
generate_plots = false;

EEG = simulate_phaseLag_data(sim_freq, sim_phaselag, dipole_1_loc, dipole_2_loc, ...
    activation_win, 'noise_level', 10, 'signal_amp', [10, 10], 'show_plots', false);

%% Preprocessing: Filter data around target frequency
filter_range = [12 18]; % Hz around 15 Hz
nyquist = EEG.srate/2;
frange = filter_range / nyquist;
transw = 0.1;
order = round(3 * EEG.srate / filter_range(1));

% Design filter
shape = [0 0 1 1 0 0];
frex = [0 frange(1)-frange(1)*transw frange frange(2)+frange(2)*transw 1];
filtkern = firls(order, frex, shape);

% Apply filter to all trials
EEG.fdata = zeros(EEG.nbchan, EEG.pnts, EEG.trials);
for triali = 1:EEG.trials
    for chani = 1:EEG.nbchan
        EEG.fdata(chani,:,triali) = filtfilt(filtkern, 1, EEG.data(chani,:,triali));
    end
end

fprintf('Data filtered between %d-%d Hz\n', filter_range);

%% Compute covariance matrices for beamforming
active_time = dsearchn(EEG.times', activation_win');
baseline_time = dsearchn(EEG.times', [-0.5 -0.1]');

% Active period covariance
S_active = zeros(EEG.nbchan);
for triali = 1:EEG.trials
    tmp = EEG.fdata(:, active_time(1):active_time(2), triali);
    tmp = bsxfun(@minus, tmp, mean(tmp,2));
    S_active = S_active + (tmp * tmp') / size(tmp,2);
end
S_active = S_active / EEG.trials;

% Full data covariance (for regularization)
S_full = zeros(EEG.nbchan);
for triali = 1:EEG.trials
    tmp = EEG.fdata(:, :, triali);
    tmp = bsxfun(@minus, tmp, mean(tmp,2));
    S_full = S_full + (tmp * tmp') / size(tmp,2);
end
S_full = S_full / EEG.trials;

%% Standard LCMV Beamformer with improved regularization
n_dipoles = size(lf.Gain, 3);
neural_ai = zeros(n_dipoles, 1);
source_power = zeros(n_dipoles, 1);

% Regularization parameter
reg_param = 0.05 * trace(S_full) / EEG.nbchan;
S_inv = inv(S_full + reg_param * eye(EEG.nbchan));

% Compute beamformer weights and neural activity index
for di = 1:n_dipoles
    lf_dip = squeeze(lf.Gain(:, 1, di));
    
    % LCMV spatial filter
    denom = lf_dip' * S_inv * lf_dip;
    if denom > eps
        w = (S_inv * lf_dip) / denom;
        
        % Source power during active period
        source_power(di) = real(w' * S_active * w);
        
        % Neural activity index (pseudo-Z score)
        neural_ai(di) = source_power(di);
    end
end

% Normalize
neural_ai = neural_ai / max(neural_ai);

%% Alternative: DICS beamformer for better handling of correlated sources
% Cross-spectral density at target frequency
csd_active = zeros(EEG.nbchan);
csd_baseline = zeros(EEG.nbchan);

% Compute CSD using Fourier transform
nfft = 2^nextpow2(diff(active_time)+1);
freq_idx = round(sim_freq * nfft / EEG.srate) + 1;

for triali = 1:EEG.trials
    % Active period
    tmp = EEG.fdata(:, active_time(1):active_time(2), triali);
    fft_data = fft(tmp, nfft, 2);
    csd_active = csd_active + (fft_data(:,freq_idx) * fft_data(:,freq_idx)') / EEG.trials;
    
    % Baseline
    tmp = EEG.fdata(:, baseline_time(1):baseline_time(2), triali);
    fft_data = fft(tmp, nfft, 2);
    csd_baseline = csd_baseline + (fft_data(:,freq_idx) * fft_data(:,freq_idx)') / EEG.trials;
end

% DICS spatial filters
dics_power = zeros(n_dipoles, 1);
reg_param_dics = 0.1 * trace(real(csd_active)) / EEG.nbchan;
csd_inv = inv(csd_active + reg_param_dics * eye(EEG.nbchan));

for di = 1:n_dipoles
    lf_dip = squeeze(lf.Gain(:, 1, di));
    
    % DICS filter
    denom = real(lf_dip' * csd_inv * lf_dip);
    if denom > eps
        w = (csd_inv * lf_dip) / denom;
        dics_power(di) = real(w' * csd_active * w);
    end
end

dics_power = dics_power / max(dics_power);

%% Combine LCMV and DICS results
combined_metric = 0.5 * neural_ai + 0.5 * dics_power;

%% Find peaks with spatial constraint
[sorted_metric, sort_idx] = sort(combined_metric, 'descend');

% Safety check
if length(sort_idx) < 2
    error('Not enough dipoles detected. Check your data and parameters.');
end

% Display results
fprintf('\nTop 10 combined metric peaks:\n');
for i = 1:min(10, length(sort_idx))
    dist1 = sqrt(sum((lf.GridLoc(sort_idx(i),:) - lf.GridLoc(dipole_1_loc,:)).^2));
    dist2 = sqrt(sum((lf.GridLoc(sort_idx(i),:) - lf.GridLoc(dipole_2_loc,:)).^2));
    fprintf('  Peak %d: Dipole %d, Metric = %.3f, Dist to sim: %.1f/%.1f mm\n', ...
        i, sort_idx(i), sorted_metric(i), dist1, dist2);
end

% Find reconstructed dipoles
% First peak
recon_dip1 = sort_idx(1);

% Second peak: must be spatially separated
min_separation = 20; % mm
recon_dip2 = []; % Initialize

for i = 2:length(sort_idx)
    dist = sqrt(sum((lf.GridLoc(sort_idx(i),:) - lf.GridLoc(recon_dip1,:)).^2));
    if dist > min_separation
        recon_dip2 = sort_idx(i);
        break;
    end
end

% Fallback if no dipole meets separation criterion
if isempty(recon_dip2)
    fprintf('Warning: No dipole found with minimum separation. Using dipole #%d\n', sort_idx(2));
    recon_dip2 = sort_idx(2);
end

%% Evaluate reconstruction accuracy
dist_error1 = min(sqrt(sum((lf.GridLoc(recon_dip1,:) - lf.GridLoc(dipole_1_loc,:)).^2)), ...
                  sqrt(sum((lf.GridLoc(recon_dip1,:) - lf.GridLoc(dipole_2_loc,:)).^2)));
dist_error2 = min(sqrt(sum((lf.GridLoc(recon_dip2,:) - lf.GridLoc(dipole_1_loc,:)).^2)), ...
                  sqrt(sum((lf.GridLoc(recon_dip2,:) - lf.GridLoc(dipole_2_loc,:)).^2)));

fprintf('\n=== RECONSTRUCTION RESULTS ===\n');
fprintf('Original dipoles: %d, %d\n', dipole_1_loc, dipole_2_loc);
fprintf('Reconstructed dipoles: %d, %d\n', recon_dip1, recon_dip2);
fprintf('Localization errors: %.1f mm, %.1f mm\n', dist_error1, dist_error2);
fprintf('Combined metric values: %.3f, %.3f\n', ...
    combined_metric(recon_dip1), combined_metric(recon_dip2));

%% Visualize results
figure('Position', [100 100 1200 400]);

% Original dipole projections
subplot(141)
topoplotIndie(-lf.Gain(:,1,dipole_1_loc), EEG.chanlocs, 'numcontour', 0);
title(sprintf('Original Dipole 1 (#%d)', dipole_1_loc))

subplot(142)
topoplotIndie(-lf.Gain(:,1,dipole_2_loc), EEG.chanlocs, 'numcontour', 0);
title(sprintf('Original Dipole 2 (#%d)', dipole_2_loc))

% Reconstructed dipole projections
subplot(143)
topoplotIndie(-lf.Gain(:,1,recon_dip1), EEG.chanlocs, 'numcontour', 0);
title(sprintf('Reconstructed 1 (#%d)\nError: %.1f mm', recon_dip1, dist_error1))

subplot(144)
topoplotIndie(-lf.Gain(:,1,recon_dip2), EEG.chanlocs, 'numcontour', 0);
title(sprintf('Reconstructed 2 (#%d)\nError: %.1f mm', recon_dip2, dist_error2))

%% Reconstruct source time series
% Use better of LCMV or DICS weights
if mean([neural_ai(recon_dip1), neural_ai(recon_dip2)]) > ...
   mean([dics_power(recon_dip1), dics_power(recon_dip2)])
    fprintf('\nUsing LCMV weights for source reconstruction\n');
    use_lcmv = true;
else
    fprintf('\nUsing DICS-based approach for source reconstruction\n');
    use_lcmv = false;
end

source_ts1 = zeros(EEG.pnts, EEG.trials);
source_ts2 = zeros(EEG.pnts, EEG.trials);

% Get leadfield vectors
lf1 = squeeze(lf.Gain(:, 1, recon_dip1));
lf2 = squeeze(lf.Gain(:, 1, recon_dip2));

if use_lcmv
    % LCMV spatial filters
    w1 = (S_inv * lf1) / (lf1' * S_inv * lf1);
    w2 = (S_inv * lf2) / (lf2' * S_inv * lf2);
else
    % Use simple projection for frequency-domain approach
    % Normalize leadfield vectors
    w1 = lf1 / norm(lf1);
    w2 = lf2 / norm(lf2);
end

% Apply spatial filters
for triali = 1:EEG.trials
    source_ts1(:, triali) = w1' * EEG.fdata(:, :, triali);
    source_ts2(:, triali) = w2' * EEG.fdata(:, :, triali);
end

%% Compute phase synchronization
% Extract phase using Hilbert transform
phase1 = angle(hilbert(source_ts1));
phase2 = angle(hilbert(source_ts2));

% Compute ISPC over time
time_window = 100; % samples
ispc_time = zeros(1, EEG.pnts);

for t = 1:EEG.pnts
    t_start = max(1, t - time_window/2);
    t_end = min(EEG.pnts, t + time_window/2);
    
    phase_diff = phase1(t_start:t_end, :) - phase2(t_start:t_end, :);
    ispc_time(t) = abs(mean(exp(1i * phase_diff(:))));
end

% Plot synchronization
figure('Position', [100 100 800 600]);

subplot(311)
plot(EEG.times, mean(source_ts1, 2), 'b', 'LineWidth', 1.5);
hold on
plot(EEG.times, mean(source_ts2, 2), 'r', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Amplitude');
title('Reconstructed Source Time Series');
legend('Source 1', 'Source 2');
xlim([-0.5 2.5]);

subplot(312)
plot(EEG.times, ispc_time, 'k', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('ISPC');
title('Phase Synchronization Over Time');
xlim([-0.5 2.5]); ylim([0 1]);

subplot(313)
imagesc(EEG.times, 1:2, [combined_metric(recon_dip1), combined_metric(recon_dip2)]' * ones(1, EEG.pnts));
xlabel('Time (s)'); ylabel('Source');
title('Source Activity Strength');
colorbar; xlim([-0.5 2.5]);