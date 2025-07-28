%% GED Tutorial with Simulated EEG Data
clear; close all; clc

%% 1. Setup EEG structure and parameters
load emptyEEG  % Contains leadfield and channel locations

EEG.srate = 500;
EEG.pnts = 1000;
EEG.times = (0:EEG.pnts-1)/EEG.srate;
EEG.trials = 100;
EEG.nbchan = 64;
EEG.data = zeros(EEG.nbchan, EEG.pnts, EEG.trials);

% Stimulus onset time
stim_onset = 1; % seconds
stim_idx = dsearchn(EEG.times', stim_onset);

% SNR parameters
signal_amplitude = 2;
noise_amplitude = 2;
snr_ratio = signal_amplitude / noise_amplitude;

%% 2. Select dipole locations
dipole1 = 109; 
dipole2 = 201;  

% Visualize dipole projections (ground truth)
figure(1), clf
subplot(221)
topoplotIndie(-lf.Gain(:,1,dipole1), EEG.chanlocs, 'numcontour', 0);
title('Ground truth: Dipole 1 (10 Hz)')

subplot(222)
topoplotIndie(-lf.Gain(:,1,dipole2), EEG.chanlocs, 'numcontour', 0);
title('Ground truth: Dipole 2 (15 Hz)')

%% 3. Simulate dipole activity
freq1 = 10; % Hz - active after stimulus
freq2 = 15; % Hz - active after stimulus

for trial = 1:EEG.trials
    % Background noise for all dipoles
    dipole_activity = randn(size(lf.Gain,3), EEG.pnts) * 0.1;
    
    % Dipole 1: active after stimulus (10 Hz)
    pre_stim = zeros(1, stim_idx);
    post_stim = signal_amplitude * sin(2*pi*freq1*EEG.times(stim_idx+1:end) + rand*2*pi);
    dipole_activity(dipole1,:) = [pre_stim, post_stim];
    
    % Dipole 2: active after stimulus (15 Hz)
    pre_stim = zeros(1, stim_idx);
    post_stim = signal_amplitude * sin(2*pi*freq2*EEG.times(stim_idx+1:end) + rand*2*pi);
    dipole_activity(dipole2,:) = [pre_stim, post_stim];
    
    % Project to scalp
    EEG.data(:,:,trial) = squeeze(lf.Gain(:,1,:)) * dipole_activity;
end

% Add sensor noise
EEG.data = EEG.data + randn(size(EEG.data)) * noise_amplitude;

%% 4. Visualize simulated data
figure(2), clf
subplot(211)
plot(EEG.times, mean(EEG.data(31,:,:), 3))
hold on
plot([stim_onset stim_onset], ylim, 'r--', 'LineWidth', 2)
title('Channel 31 - Average ERP')
xlabel('Time (s)'), ylabel('Amplitude')

subplot(212)
imagesc(EEG.times, 1:EEG.trials, squeeze(EEG.data(31,:,:))')
hold on
plot([stim_onset stim_onset], ylim, 'r--', 'LineWidth', 2)
title('Channel 31 - Single trials')
xlabel('Time (s)'), ylabel('Trial')

%% 5. Define time windows for GED
pre_window = [0.2 0.8];     % pre-stimulus
post_window = [1.2 1.8];    % post-stimulus

pre_idx = dsearchn(EEG.times', pre_window');
post_idx = dsearchn(EEG.times', post_window');

% Highlight windows
figure(3), clf
plot(EEG.times, mean(mean(EEG.data,3),1))
hold on
% Get the y-axis limits *after* plotting the data
yl = ylim; 
% Create the patches using the retrieved limits
patch([pre_window fliplr(pre_window)], [yl(1) yl(1) yl(2) yl(2)], 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
patch([post_window fliplr(post_window)], [yl(1) yl(1) yl(2) yl(2)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot([stim_onset stim_onset], ylim, 'k--', 'LineWidth', 2)
legend('Global signal', 'Pre-stim', 'Post-stim', 'Stimulus')
xlabel('Time (s)'), ylabel('Amplitude')
title('Time windows for covariance matrices')

%% 6. Compute covariance matrices
% Pre-stimulus covariance
data_pre = EEG.data(:, pre_idx(1):pre_idx(2), :);  % extract pre-stim data
data_pre = reshape(data_pre, EEG.nbchan, []);      % concatenate trials to make it channel x time/trial
data_pre = data_pre - mean(data_pre, 2);           % mean center
cov_pre = (data_pre * data_pre') / (size(data_pre,2) - 1);  % compute covariance

% Post-stimulus covariance
data_post = EEG.data(:, post_idx(1):post_idx(2), :);  % extract post-stim data
data_post = reshape(data_post, EEG.nbchan, []);       % concatenate trials
data_post = data_post - mean(data_post, 2);           % remove mean
cov_post = (data_post * data_post') / (size(data_post,2) - 1);  % compute covariance

% Visualize covariance matrices
figure(4), clf
clim_val = max(abs([cov_pre(:); cov_post(:)])) * 0.8;

subplot(121)
imagesc(cov_pre)
axis square
colorbar
caxis([-clim_val clim_val])
title('Pre-stimulus covariance')
xlabel('Channel'), ylabel('Channel')

subplot(122)
imagesc(cov_post)
axis square
colorbar
caxis([-clim_val clim_val])
title('Post-stimulus covariance')
xlabel('Channel'), ylabel('Channel')

%% 7. Perform GED
[V, D] = eig(cov_post, cov_pre);  % solve generalized eigenvalue problem
[eigenvalues, sort_idx] = sort(diag(D), 'descend');  % sort by eigenvalue
V = V(:, sort_idx);  % reorder eigenvectors

% Plot eigenspectrum
figure(5), clf
subplot(211)
plot(eigenvalues(1:20), 'o-', 'LineWidth', 2, 'MarkerSize', 8)
xlabel('Component'), ylabel('Eigenvalue')
title('GED eigenspectrum')
grid on

%% 8. Extract and visualize components
n_comps = 2;
component_maps = zeros(n_comps, EEG.nbchan);
component_ts = zeros(n_comps, EEG.pnts, EEG.trials);

for comp = 1:n_comps
    % Component map (forward model)
    component_maps(comp,:) = V(:,comp)' * cov_post;  % compute spatial filter projection
    
    % Fix sign based on maximum
    [~, max_idx] = max(abs(component_maps(comp,:)));  % find strongest channel
    component_maps(comp,:) = component_maps(comp,:) * sign(component_maps(comp,max_idx));  % ensure positive peak
    
    % Component time series
    comp_data = V(:,comp)' * reshape(EEG.data, EEG.nbchan, []);  % apply spatial filter
    component_ts(comp,:,:) = reshape(comp_data, EEG.pnts, EEG.trials);  % reshape to trials
end

% Plot component maps (compare with ground truth)
figure(1)
subplot(223)
topoplotIndie(component_maps(1,:), EEG.chanlocs, 'numcontour', 0);
title('GED: Component 1')

subplot(224)
topoplotIndie(component_maps(2,:), EEG.chanlocs, 'numcontour', 0);
title('GED: Component 2')

%% 9. Component time series analysis
figure(6), clf
colors = {'b', 'r'};

for comp = 1:n_comps
    subplot(2,1,comp)
    
    % Plot average time series
    avg_ts = mean(component_ts(comp,:,:), 3);
    plot(EEG.times, avg_ts, colors{comp}, 'LineWidth', 2)
    hold on
    
    % Mark stimulus onset
    plot([stim_onset stim_onset], ylim, 'k--', 'LineWidth', 2)

    % Highlight analysis windows
    yl = ylim;
    % Create the patches using the retrieved limits
    patch([pre_window fliplr(pre_window)], [yl(1) yl(1) yl(2) yl(2)], 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
    patch([post_window fliplr(post_window)], [yl(1) yl(1) yl(2) yl(2)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none')

    title(['Component ' num2str(comp) ' time series'])
    xlabel('Time (s)'), ylabel('Amplitude')
    grid on
end

%% 10. Power spectral analysis of components
figure(7), clf
nfft = 2^nextpow2(EEG.pnts);
hz = linspace(0, EEG.srate/2, floor(nfft/2)+1);

for comp = 1:n_comps
    % Pre-stimulus power
    pre_data = squeeze(component_ts(comp, pre_idx(1):pre_idx(2), :));
    pre_fft = fft(pre_data, nfft, 1);
    pre_pow = mean(abs(pre_fft(1:length(hz),:)).^2, 2);
    
    % Post-stimulus power
    post_data = squeeze(component_ts(comp, post_idx(1):post_idx(2), :));
    post_fft = fft(post_data, nfft, 1);
    post_pow = mean(abs(post_fft(1:length(hz),:)).^2, 2);
    
    subplot(1,2,comp)
    plot(hz, pre_pow, 'g', 'LineWidth', 2)
    hold on
    plot(hz, post_pow, 'r', 'LineWidth', 2)
    xlim([0 40])
    xlabel('Frequency (Hz)'), ylabel('Power')
    title(['Component ' num2str(comp) ' spectrum'])
    legend('Pre-stim', 'Post-stim')
    grid on
end

%% Summary
fprintf('\nGED Analysis Summary:\n')
fprintf('-------------------\n')
fprintf('SNR ratio: %.2f\n', snr_ratio)
fprintf('Component 1: Eigenvalue = %.2f\n', eigenvalues(1))
fprintf('Component 2: Eigenvalue = %.2f\n', eigenvalues(2))
fprintf('Ratio (λ1/λ2) = %.2f\n', eigenvalues(1)/eigenvalues(2))