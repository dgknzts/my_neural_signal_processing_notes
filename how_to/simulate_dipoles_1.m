%%
%     TUTORIAL: How to Simulate Dipoles in EEG - Complete Beginner's Guide
%    
%    This tutorial will teach you:
%    1. What are dipoles in brain activity?
%    2. How to create simple dipole signals
%    3. How to project dipoles to EEG electrodes
%    4. How to simulate realistic EEG data from dipoles
%
%%

%% PART 1: UNDERSTANDING DIPOLES
% 
% What is a dipole?
% ================
% In EEG, a "dipole" represents a small group of neurons firing together.
% Think of it as a tiny electrical source inside your brain.
% 
% Key concepts:
% - Location: WHERE in the brain the activity happens (x,y,z coordinates)
% - Orientation: Which DIRECTION the electrical activity points
% - Strength: HOW STRONG the activity is over time
% - Forward Model: How brain activity projects to scalp electrodes

fprintf('=== DIPOLE SIMULATION TUTORIAL ===\n\n');

%% PART 2: BASIC SETUP
% Let's start with the simplest possible example

% First, let's create some basic parameters
srate = 500;  % sampling rate (500 Hz)
time_duration = 2;  % 2 seconds of data
times = -1:1/srate:1;  % time vector from -1 to +1 seconds
n_timepoints = length(times);

fprintf('Created time vector with %d points at %d Hz sampling rate\n', n_timepoints, srate);

%% PART 3: SIMULATE A SINGLE DIPOLE SIGNAL
% Let's create the simplest dipole: a sine wave burst

% Parameters for our dipole signal
dipole_frequency = 10;  % 10 Hz oscillation (alpha rhythm)
burst_center = 0.2;       % burst happens at time 0
burst_width = 0.3;      % burst lasts 0.3 seconds

% Create a Gaussian envelope (the "burst" shape)
gaussian_envelope = exp(-4*log(2)*(times - burst_center).^2 / burst_width^2);

% Create the oscillating signal
sine_wave = sin(2*pi * dipole_frequency * times);

% Combine them: oscillation × envelope = burst of activity
dipole_signal = sine_wave .* gaussian_envelope;

% Plot to visualize
figure(1); clf;
subplot(3,1,1)
plot(times, gaussian_envelope, 'b-', 'LineWidth', 2);
title('Step 1: Gaussian Envelope (burst shape)');
xlabel('Time (s)'); ylabel('Amplitude');
grid on;

subplot(3,1,2)
plot(times, sine_wave, 'r-');
title('Step 2: Sine Wave (10 Hz oscillation)');
xlabel('Time (s)'); ylabel('Amplitude');
grid on;

subplot(3,1,3)
plot(times, dipole_signal, 'k-', 'LineWidth', 2);
title('Step 3: Final Dipole Signal = Envelope × Sine Wave');
xlabel('Time (s)'); ylabel('Amplitude');
grid on;

fprintf('Created a dipole signal: 10 Hz burst at time 0\n');

%% PART 4: FROM BRAIN TO SCALP - THE FORWARD MODEL CONCEPT
%
% The brain signal we just created needs to travel from inside the brain
% to the scalp electrodes. This is what the "forward model" does.
%
% Think of it like this:
% Brain activity (dipole) → travels through skull/tissue → detected at scalp (EEG)

% Let's simulate a simple 3-electrode setup
n_electrodes = 3;
electrode_names = {'Fz', 'Cz', 'Pz'};  % Front, Center, Back of head

% Simple forward model: how much each electrode sees our dipole
% (In reality, this comes from solving physics equations)
forward_weights = [0.8, 1.0, 0.3];  % Center electrode sees it most

% Project our dipole signal to each electrode
eeg_data = zeros(n_electrodes, n_timepoints);
for elec = 1:n_electrodes
    eeg_data(elec, :) = forward_weights(elec) * dipole_signal;
end

% Plot the EEG data
figure(2); clf;
for elec = 1:n_electrodes
    subplot(n_electrodes, 1, elec);
    plot(times, eeg_data(elec, :), 'LineWidth', 2);
    title(sprintf('EEG at electrode %s (weight = %.1f)', electrode_names{elec}, forward_weights(elec)));
    xlabel('Time (s)'); ylabel('Amplitude (µV)');
    grid on;
end

fprintf('Projected dipole to %d electrodes using simple forward model\n', n_electrodes);

%% PART 5: ADDING NOISE (REALISTIC EEG)
%
% Real EEG is never perfectly clean. Let's add some realistic noise.

% Types of noise in EEG:
% 1. Background brain activity (pink noise)
% 2. Muscle artifacts (high frequency)
% 3. Electrical interference (50/60 Hz)

noise_level = 2;  % Adjust this to make noise stronger/weaker

eeg_data_noisy = zeros(size(eeg_data));
for elec = 1:n_electrodes
    % Pink noise (typical background EEG)
    pink_noise = noise_level * cumsum(randn(1, n_timepoints)) / sqrt(n_timepoints);
    
    % Add to clean signal
    eeg_data_noisy(elec, :) = eeg_data(elec, :) + pink_noise;
end

% Compare clean vs noisy
figure(3); clf;
subplot(2,1,1)
plot(times, eeg_data(2, :), 'b-', 'LineWidth', 2);
title('Clean EEG Signal (Center electrode)');
xlabel('Time (s)'); ylabel('Amplitude (µV)');
grid on;

subplot(2,1,2)
plot(times, eeg_data_noisy(2, :), 'r-', 'LineWidth', 1);
hold on;
plot(times, eeg_data(2, :), 'b-', 'LineWidth', 2);
title('Noisy EEG Signal (red) vs Clean Signal (blue)');
xlabel('Time (s)'); ylabel('Amplitude (µV)');
legend('Noisy', 'Clean', 'Location', 'best');
grid on;

fprintf('Added realistic noise to EEG signals\n');

%% PART 6: MULTIPLE DIPOLES
%
% Real brains have many active regions. Let's simulate two dipoles.

% Dipole 1: 10 Hz activity in "frontal" region
dipole1_freq = 10;
dipole1_envelope = exp(-4*log(2)*(times - 0.2).^2 / 0.3^2);
dipole1_signal = sin(2*pi * dipole1_freq * times) .* dipole1_envelope;
dipole1_weights = [1.0, 0.5, 0.1];  % Strong in front

% Dipole 2: 20 Hz activity in "parietal" region  
dipole2_freq = 20;
dipole2_envelope = exp(-4*log(2)*(times + 0.2).^2 / 0.2^2);
dipole2_signal = sin(2*pi * dipole2_freq * times) .* dipole2_envelope;
dipole2_weights = [0.1, 0.5, 1.0];  % Strong in back

% Combine both dipoles at each electrode
eeg_multi = zeros(n_electrodes, n_timepoints);
for elec = 1:n_electrodes
    eeg_multi(elec, :) = dipole1_weights(elec) * dipole1_signal + ...
                         dipole2_weights(elec) * dipole2_signal;
    
    % Add noise
    noise = 0.1 * cumsum(randn(1, n_timepoints)) / sqrt(n_timepoints);
    eeg_multi(elec, :) = eeg_multi(elec, :) + noise;
end

% Plot the multi-dipole result
figure(4); clf;
for elec = 1:n_electrodes
    subplot(n_electrodes, 1, elec);
    plot(times, eeg_multi(elec, :), 'LineWidth', 1.5);
    title(sprintf('Multi-dipole EEG at %s', electrode_names{elec}));
    xlabel('Time (s)'); ylabel('Amplitude (µV)');
    grid on;
end

fprintf('Simulated EEG from two dipoles with different frequencies and locations\n');

%% PART 7: MULTIPLE TRIALS (EXPERIMENTAL DESIGN)
%
% In real experiments, we record many trials and average them.

n_trials = 50;
trial_data = zeros(n_electrodes, n_timepoints, n_trials);

fprintf('Simulating %d trials...\n', n_trials);

for trial = 1:n_trials
    % Add some trial-to-trial variability
    
    % Dipole 1: frequency varies slightly across trials
    freq_jitter = dipole1_freq + randn*10;  % 10 ± 0.5 Hz
    trial_dipole1 = sin(2*pi * freq_jitter * times) .* dipole1_envelope;
    
    % Dipole 2: timing varies slightly
    time_jitter = 0.05 * randn;  % ± 50ms timing jitter
    shifted_envelope = exp(-4*log(2)*(times + 0.2 - time_jitter).^2 / 0.2^2);
    trial_dipole2 = sin(2*pi * dipole2_freq * times) .* shifted_envelope;
    
    % Project to electrodes and add noise
    for elec = 1:n_electrodes
        signal = dipole1_weights(elec) * trial_dipole1 + ...
                dipole2_weights(elec) * trial_dipole2;
        noise = 0.15 * cumsum(randn(1, n_timepoints)) / sqrt(n_timepoints);
        trial_data(elec, :, trial) = signal + noise;
    end
end

% Calculate average across trials
average_data = mean(trial_data, 3);

% Plot single trial vs average
figure(5); clf;
subplot(2,1,1)
plot(times, squeeze(trial_data(2, :, 1)), 'Color', [0.7 0.7 0.7]);
hold on;
plot(times, squeeze(trial_data(2, :, 2)), 'Color', [0.7 0.7 0.7]);
plot(times, squeeze(trial_data(2, :, 3)), 'Color', [0.7 0.7 0.7]);
title('Single Trials (gray) - Center Electrode');
xlabel('Time (s)'); ylabel('Amplitude (µV)');
grid on;

subplot(2,1,2)
plot(times, average_data(2, :), 'b-', 'LineWidth', 3);
title(sprintf('Average of %d Trials - Noise Reduced!', n_trials));
xlabel('Time (s)'); ylabel('Amplitude (µV)');
grid on;

fprintf('Trial averaging reduces noise and reveals true signal structure\n');


%% PARAMETER EXPLORATION
%
% Try changing these parameters to see their effects:

% Experiment with these:
% - dipole_frequency: Change from 10 to other values (1-100 Hz)
% - burst_width: Make bursts shorter (0.1) or longer (0.5)
% - noise_level: Less noise (0.05) or more noise (0.3)
% - forward_weights: Different spatial patterns
% - n_trials: See how averaging improves with more trials
