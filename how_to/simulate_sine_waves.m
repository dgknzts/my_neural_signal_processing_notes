%% Sine wave
% formula of a basic sine wave is:
% amplitude * sin(2*pi* frequency* time + phase )
% amplitude : representing the y axis -> kinda strength of the wave
% frequency : number of cycles within a time unit (second)
% time : length of the signal/wave representing sampling rate and length
% phase : location of the wave global shift of the wave to the left or right

% note : if phase is zero (0) its a sine wave
% if phase is pi/2 cosine wave. there is exactly 0 correlation
% between them.

% intuition about sine waves
% define some variables
freq  = 5;    % frequency in Hz
srate = 500; % sampling rate in Hz
time  = -1:1/srate:1; % time vector in seconds
ampl  = 1;
phas  = 0;%pi/3;

% now create the sinewave
sinewave = ampl*sin(2*pi*freq*time + phas );

% now plot it!
figure(1), clf
plot(time,sinewave)

% optional prettification of the plot
set(gca,'xlim',[min(time-0.1) max(time) + 0.1],'ylim',[-1.1 1.1]*ampl); 
xlabel('Time (ms)')
title('My first sine wave plot! Mom will be so proud!')

%% The sum of sine waves can appear like a complicated time series
% Create a few sine waves and sum
% List some frequencies
% Parameters
frex = [4.8, 6, 7.5];  % Frequencies
amplit = [5, 5, 5];    % Amplitudes
phases = [2*pi, 2*pi, 2*pi];  % Phases
srate = 1000;  % Sampling rate (Hz)
time = 0:1/srate:1;  % Time vector

% Generate sine waves
sine_waves = zeros(length(frex), length(time));
for fi = 1:length(frex)
    sine_waves(fi,:) = amplit(fi) * sin(2*pi*frex(fi)*time + phases(fi));
    sine_waves(fi,:) = (sine_waves(fi,:) - min(sine_waves(fi,:))) / (max(sine_waves(fi,:)) - min(sine_waves(fi,:))) * 0.5 + 0.5; % Normalize
end

% Plot
figure, clf
colors = ["#F4A261", "#1C646D", "#38184C"];  % Custom colors
arc_labels = ["Inner Arc", "Middle Arc", "Outer Arc"];  % Labels

% Plot individual sine waves
for fi = 1:length(frex)
    subplot(length(frex)+1, 1, fi)
    plot(time, sine_waves(fi,:), 'Color', colors(fi), 'LineWidth', 1.5)
    ylim([0.5, 1]), yticks([0.5, 0.75, 1])
    title(sprintf('%s: %.1f Hz', arc_labels(fi), frex(fi)))
    ylabel('Contrast')
end

% Plot the sum of sine waves
subplot(length(frex)+1, 1, length(frex)+1)
plot(time, sum(sine_waves), 'k', 'LineWidth', 1.5)
title('Sum of Sine Waves')
xlabel('Time (s)'), ylabel('Amplitude')
