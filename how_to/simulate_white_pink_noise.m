%%%% How to simulate neural noise %%%%
%% White Noise: Mostly effect the time domain signal. 
% Time and frequency domain signal looks similar.
% Normally distributed white noise.

% simulation details
srate  = 200; % sampling rate in Hz
time   = -1:1/srate:2;
pnts   = length(time);
signal = 10*sin(2*pi*time*20); % signal as a pure sine wave
% frequencies for the power spectrum
hz = linspace(0,srate/2,floor(length(time)/2)-1);
% time/2 nyquist theorem


% Parameters:
stretch = 10; % wideness of the gaussian
shift   = -5; % shift from zero amplitude

rng(3); % works like a SEED

% generate random data
white_noise = stretch*randn(size(time)) + shift;

%% Pink (1/f or fractal) Noise: Mostly effect the frequency domain signal
% This is more similar to real neural data compared to white.
rng(12)

% Parameters:
% generate 1/f amplitude spectrum
ed = 50; % exponential decay parameter. in other words the slope of the distribution
as = 3*rand(1,floor(pnts/2)-1) .* exp(-(1:floor(pnts/2)-1)/ed); 
%random amplitude values multiplied by exponentially decreasing decay parameter generator
%noise is decreasing for higher frequencies. we have to apply this noise in
%frequency domain than convert it back to the time domain because pink
%noise is related to the frequency domain activity

%trying out the fixed amplitude to see the decrease rate clean
%as = exp(-(1:floor(pnts/2)-1)/ed); 

as = [as(1) as 0 0 as(:,end:-1:1)]; % mirroring just like in the fft. neg and pos values
% more clear after fft 

% Fourier coefficients
fc = as .* exp(1i*2*pi*rand(size(as))); %second equation is adding random phase values.

% inverse Fourier transform to create the noise
pink_noise = real(ifft(fc)) * pnts; %inverse fft normalize the signal by dividing it with n automatically
% so we multipyle it back to the real signal.

%% Plotting the noise
% Play with the upper parameters to see the difference
% Compute frequency domain representations
signal_fft = 2*abs(fft(signal) / pnts);
white_noise_fft = 2*abs(fft(signal + white_noise)/ pnts);
pink_noise_fft = 2*abs(fft(signal + pink_noise)/ pnts);
both_noise_fft = 2*abs(fft(signal + white_noise + pink_noise)/ pnts);

% Define signal variations and their properties
signals = {signal, signal_fft(1:length(hz)), ...
           signal + white_noise, white_noise_fft(1:length(hz)), ...
           signal + pink_noise, pink_noise_fft(1:length(hz)), ...
           signal + white_noise + pink_noise, both_noise_fft(1:length(hz))};

titles = {'Pure Signal (Time Domain)', 'Pure Signal (Frequency Domain)', ...
          'Signal with White Noise (Time Domain)', 'Signal with White Noise (Frequency Domain)', ...
          'Signal with Pink Noise (Time Domain)', 'Signal with Pink Noise (Frequency Domain)', ...
          'Signal with Both Noises (Time Domain)', 'Signal with Both Noises (Frequency Domain)'};

colors = {'k', 'k', 'r', 'r', 'g', 'g', 'b', 'b'}; % Colors for each plot
x_labels = {'Time (s)', 'Frequency (Hz)'}; % X-axis labels

% Define X-axis limits based on parameters
time_xlim = [min(time), max(time)];
freq_xlim = [0, srate/2];

% Compute global y-axis limits dynamically
time_ylim = [min([signals{1:2:7}]) - 2, max([signals{1:2:7}]) + 2]; % For time-domain plots
freq_ylim = [0, max([signals{2:2:8}]) + 0.1 * max([signals{2:2:8}])]; % For frequency-domain plots

% Create figure and plot
figure;
for i = 1:8
    subplot(4,2,i);
    if mod(i,2) == 1
        plot(time, signals{i}, colors{i}); % Time domain
        xlabel(x_labels{1});
        xlim(time_xlim);
        ylim(time_ylim);
    else
        plot(hz, signals{i}, colors{i}); % Frequency domain
        xlabel(x_labels{2});
        xlim(freq_xlim);
        ylim(freq_ylim);
    end
    ylabel('Amplitude');
    title(titles{i});
end

% Adjust layout
sgtitle('Neural Noise Simulation');