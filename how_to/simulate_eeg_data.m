%% simulating data is similar to a real eeg data with channel and trial dimensions
% this parameters are similary designed to the eeglab data structure
% parameters
EEG.srate  = 500; % sampling rate in Hz
EEG.pnts   = 1500;
EEG.trials = 30;
EEG.nbchan = 23;

sinefreq = 6.75; % in Hz

% time vector
EEG.times = (0:EEG.pnts-1)/EEG.srate;

% loop over channels and trials to create data
for chani=1:EEG.nbchan
    for triali=1:EEG.trials
        % data as a sine wave
        % hint: instantaneous frequency via interpolated random numbers
        freqmod = 10 + 10*interp1(rand(1,5),linspace(1,5,EEG.pnts));
        signal  = sin( 2*pi * ((EEG.times + cumsum(freqmod))/EEG.srate) );

        EEG.data(chani,:,triali) = signal; 
        % you can adapt the signal here
        % could combine it with noise or some transient activity etc from
        % other simulations scripts
    end
end


% plot time and freq domain
addpath(genpath('src'))
plot_simEEG(EEG,2,2) % a basic function to plot different domains of the data.
% data , channel to plot, figure

%% another example with some pink noise and transient frequency
% amount of noise
noiseamp = .1;


peaktime = 1; % seconds
width    = .12;
sinefreq = 7; % for sine wave

% create Gaussian taper
gaus = exp( -(EEG.times-peaktime).^2 / (2*width^2) );


% loop over channels and trials
for chani=1:EEG.nbchan
    for triali=1:EEG.trials
        
        % trial-unique sine wave
        cosw = cos(2*pi*sinefreq*EEG.times + 2*pi*rand);
        
        %%% 1/f noise
        ed = 50; % exponential decay parameter
        as = rand(1,floor(EEG.pnts/2)-1) .* exp(-(1:floor(EEG.pnts/2)-1)/ed);
        as = [as(1) as 0 as(:,end:-1:1)];
        
        % Fourier coefficients
        fc = as .* exp(1i*2*pi*rand(size(as)));
        
        % inverse Fourier transform to create the noise
        noise = real(ifft(fc)) * EEG.pnts;
        
        
        % data as signal + noise
        EEG.data(chani,:,triali) = cosw .* gaus + noiseamp*noise;
    end
end

plot_simEEG(EEG,2,7)