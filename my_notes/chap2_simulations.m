%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Chapter 2: Simulations%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('data/chap2')

%% Normally distributed noise %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% simulation details
srate  = 200; % sampling rate in Hz
time   = -1:1/srate:2;
pnts   = length(time);

% frequencies for the power spectrum
hz = linspace(0,srate/2,floor(length(time)/2)-1);
% time/2 nyquist theorem

% noise parameters
stretch = 15; % wideness
shift   = -10; % shift from zero

% (optional: fix the random-number generator state)
% rng(3); % works like a SEED

% generate random data
noise = stretch*randn(size(time)) + shift;

% plotting
figure(4), clf
subplot(211)
plot(time,noise,'k')
set(gca,'fontsize',15)
title('Normally distributed: Time domain')
xlabel('Time (s)'), ylabel('Amplitude')

subplot(223)
[y,x] = hist(noise,100);
plot(x,y,'k','linew',2)
xlabel('Values'), ylabel('N per bin')
title('Signal histogram (distribution)')
set(gca,'fontsize',15,'xlim',[min(x) max(x)])

subplot(224)
amp = abs(fft(noise)/pnts);
amp(2:end) = 2*amp(2:end);
plot(hz,amp(1:length(hz)),'k')
title('Frequency domain')
set(gca,'fontsize',15)
xlabel('Frequency (Hz)'), ylabel('Amplitude')

%% Pink noise (aka 1/f aka fractal) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulation details for this video
srate = 500; % sampling rate in Hz
time  = -1:1/srate:2;
pnts  = length(time);
hz    = linspace(0,srate/2,floor(length(time)/2)+1);

rng(12)
% generate 1/f amplitude spectrum
ed = 50; % exponential decay parameter
as = rand(1,floor(pnts/2)-1) .* exp(-(1:floor(pnts/2)-1)/ed); %random amplitude values multiplied by exponentially decreasing decay parameter generator

%trying out the fixed amplitude
%as = exp(-(1:floor(pnts/2)-1)/ed); % to see how line 71 works.

as = [as(1) as 0 0 as(:,end:-1:1)]; % mirroring just like in the fft. neg and pos values

% Fourier coefficients
fc = as .* exp(1i*2*pi*rand(size(as))); %second equation is adding random phase values.

% inverse Fourier transform to create the noise
noise = real(ifft(fc)) * pnts; %inverse fft normalize the signal by dividing it with n automatically
% so we multipyle it back to the real signal.


figure(5), clf
subplot(211)
plot(time,noise,'k')
set(gca,'fontsize',15)
title('Pink noise: Time domain')
xlabel('Time (s)'), ylabel('Amplitude')

subplot(223)
[y,x] = hist(noise,100);
plot(x,y,'k','linew',2)
xlabel('Values'), ylabel('N per bin')
title('Signal histogram (distribution)')
set(gca,'fontsize',15)

subplot(224)
amp = abs(fft(noise)/pnts);
amp(2:end) = 2*amp(2:end);
plot(hz,amp(1:length(hz)),'k')
title('Frequency domain')
set(gca,'fontsize',15)
xlabel('Frequency (Hz)'), ylabel('Amplitude')

%% Important Formulas %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Sine wave
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

% Define a sampling rate
srate = 1000;

% List some frequencies
frex = [ 4.8   6   7.5 ];
labels = ['Inner' 'Middle' 'Outer'];

% List some random amplitudes... make sure there are the same number of
% amplitudes as there are frequencies!
amplit = [5 5 5 ];

% Phases... list some random numbers between -pi and pi
phases = [2*pi 2*pi 2*pi]; 

% Define time
time = 0:1/srate:1;

% Now we loop through frequencies and create sine waves
sine_waves = zeros(length(frex),length(time));
for fi = 1:length(frex)
    sine_waves(fi,:) = amplit(fi) * sin(2*pi*time*frex(fi) + phases(fi));
    % Normalize to range [0.5, 1]
    sine_waves(fi,:) = (sine_waves(fi,:) - min(sine_waves(fi,:))) / (max(sine_waves(fi,:)) - min(sine_waves(fi,:))) * 0.5 + 0.5;
end

% Plot the sum of sine waves
figure(2), clf
plot(time, sum(sine_waves))
title('Sum of sine waves')
xlabel('Time (s)'), ylabel('Amplitude (arb. units)')

% Plot each wave separately
figure(3), clf
colors = ["#F4A261", "#1C646D", "#38184C"]; % Define custom colors for the sine waves
arc_labels = ["Inner Arc", "Middle Arc", "Outer Arc"]; % Define arc labels for each wave
for fi = 1:length(frex)
    subplot(length(frex), 1, fi)
    plot(time, sine_waves(fi,:), 'Color', colors(fi), 'LineWidth', 1.5)
    axis([time([1 end]) 0.5 1])
    ylabel('Contrast')
    yticks([0.5 0.75 1]) % Manually set x-axis ticks
    title(sprintf('%s: %.1f Hz (f%d)', arc_labels(fi), frex(fi), fi)) % Add informative title with arc labels
end
xlabel('Time (s)')


%% Gaussian

% simulation parameters
time  = -2:1/1000:2;
ptime = 1;  % peak time 
ampl  = 45; % amplitude
fwhm  = .9;

% Gaussian
gwin = ampl* exp(-(4*log(2)*(time-ptime).^2) / fwhm^2);

% empirical FWHM
gwinN   = gwin./max(gwin);
midp    = dsearchn(time',1);
pst5    = midp-1+dsearchn(gwinN(midp:end)',.5);
pre5    = dsearchn(gwinN(1:midp)',.5);
empfwhm = time(pst5) - time(pre5);


figure(3), clf, hold on
plot(time,gwin,'k','linew',2)
plot(time([pre5 pst5]),gwin([pre5 pst5]),'ro--','markerfacecolor','k')
plot(time([pre5 pre5]),[0 gwin(pre5)],'r:')
plot(time([pst5 pst5]),[0 gwin(pst5)],'r:')
title([ 'Requested FWHM: ' num2str(fwhm) 's, empirical FWHM: ' num2str(empfwhm) 's' ])
xlabel('Time (s)'), ylabel('Amplitude')

%% Euler's formula

M = 3;
k = 3*pi/2;

meik = M*exp(1i*k);

figure(4), clf
subplot(121)
polar([0 k],[0 M],'r'), hold on
polar(k, M,'ro')
title('Polar plane')


subplot(122), hold on
plot(meik,'ro')
plot(real(meik),imag(meik),'gs')
axis([-1 1 -1 1]*abs(meik))
axis square
xlabel('Real'), ylabel('Imag')
grid on
title('Cartesian (rectangular) plane')

%% chirps: simulated sine waves varying frequency acrosss time
%% simulation details

pnts  = 10000;
srate = 1024;
time  = (0:pnts-1)/srate;



% "bipolar" chirp
freqmod = linspace(5,15,pnts);

% multipolar chirp
k = 5; % poles for frequencies
freqmod = 20*interp1(rand(1,k),linspace(1,k,pnts));


% signal time series
signal  = sin( 2*pi * ((time + cumsum(freqmod))/srate) );
% Note: the code in the video has a bug in the previous line,
%   due to incorrect parentheses: srate should scale both time
%   and freqmod, as above.

% plotting

figure(1), clf

subplot(211)
plot(time,freqmod,'r','linew',3)
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('Instantaneous frequency')

subplot(212)
plot(time,signal,'k')
xlabel('Time (s)')
ylabel('Amplitude (a.u.)')
title('Signal (chirp)')

%% Non-stationary narrowband activity via filtered noise
% noise varying in time

% simulation details
pnts  = 4567;
srate =  987;

% signal parameters in Hz
peakfreq = 10; % most frequent frequency
fwhm     =  1; % bigger means frequency will vary more
% smaller for just the amplitude to vary

% frequencies
hz = linspace(0,srate,pnts);


%%% create frequency-domain Gaussian
s  = fwhm*(2*pi-1)/(4*pi); % normalized width
x  = hz-peakfreq;          % shifted frequencies
fg = exp(-.5*(x/s).^2);    % gaussian


% Fourier coefficients of random spectrum
fc = rand(1,pnts) .* exp(1i*2*pi*rand(1,pnts));

% taper with Gaussian
fc = fc .* fg;

% go back to time domain to get EEG data
signal = 2*real( ifft(fc) );


%%% plotting

figure(1), clf
subplot(211)
plot(hz,abs(fc),'k')
set(gca,'xlim',[0 peakfreq*3])
xlabel('Frequency (Hz)'), ylabel('Amplitude (a.u.)')
title('Frequency domain')


subplot(212)
plot((0:pnts-1)/srate,signal,'b')
title('Time domain')
xlabel('Time (s)'), ylabel('Amplitude')


%% Generating transient oscillations

%simulation details

pnts  = 4000;
srate = 1000;
time  = (0:pnts-1)/srate - 1;

% gaussian parameters
peaktime = 1; % seconds
fwhm     = .4;

% sine wave parameters
sinefreq = 7; % for sine wave

%% create signal

% create Gaussian taper
gaus = exp( -(4*log(2)*(time-peaktime).^2) / fwhm^2 );

% sine wave with random phase value ("non-phase-locked")
cosw = cos(2*pi*sinefreq*time + 2*pi*rand);

% signal
signal = cosw .* gaus;


% and plot
figure(1), clf
plot(time,signal,'k','linew',2)
xlabel('Time (s)'), ylabel('Amplitude')


%% 1) pure phase-locked sine wave

% parameters
EEG.srate  = 500; % sampling rate in Hz
EEG.pnts   = 1500;
EEG.trials = 30;
EEG.nbchan = 23;

sinefreq = 6.75; % in Hz

% time vector
EEG.times = (0:EEG.pnts-1)/EEG.srate;


% loop over channels and create data
for chani=1:EEG.nbchan
    for triali=1:EEG.trials
        % data as a sine wave
        EEG.data(chani,:,triali) = sin(2*pi*sinefreq*EEG.times);
    end
end


% plot an ERP from one channel
figure(1), clf
plot(EEG.times,squeeze(mean(EEG.data(10,:,:),3)),'linew',2)
xlabel('Time (s)'), ylabel('Activity')
title('ERP from channel 10')


% the function below takes at least one argument (EEG),
% and optionally a second argument (channel number),
% and optionally a third argument (figure number)
plot_simEEG(EEG,2,1)

%% 2) Non-phase-locked sine wave

% hint: copy/paste the code above but add something inside the sine 
%       function on each trial.
% parameters
EEG.srate  = 500; % sampling rate in Hz
EEG.pnts   = 1500;
EEG.trials = 30;
EEG.nbchan = 23;

sinefreq = 6.75; % in Hz

% time vector
EEG.times = (0:EEG.pnts-1)/EEG.srate;


% loop over channels and create data
for chani=1:EEG.nbchan
    for triali=1:EEG.trials
        % data as a sine wave
        EEG.data(chani,:,triali) = sin(2*pi*sinefreq*EEG.times + 2*pi*rand) ;
    end
end


% plot an ERP from one channel
figure(1), clf
plot(EEG.times,squeeze(mean(EEG.data(10,:,:),3)),'linew',2)
xlabel('Time (s)'), ylabel('Activity')
title('ERP from channel 10')


% the function below takes at least one argument (EEG),
% and optionally a second argument (channel number),
% and optionally a third argument (figure number)
plot_simEEG(EEG,2,2)

%% 3) multisine waves
rng(33)
% list of frequencies and corresponding amplitudes
frex = [ 3 5 16 ];
amps = [ 3 4 5  ];
rphase = 2*pi*rand(1, EEG.trials);

sinewave = zeros(length(frex), length(EEG.times));
% loop over channels and trials
for chani=1:EEG.nbchan
    for triali=1:EEG.trials
        for fi=1:length(frex)
            sinewave(fi,:) = amps(fi) * cos(2*pi*frex(fi)*EEG.times+ rphase(triali) );
        end
        sinewave_mix = sum(sinewave, 1);
        % data as a sine wave plus noise
        EEG.data(chani,:,triali) = sinewave_mix * randn();
    end
end
plot_simEEG(EEG,2,3)


%% 4) nonstationary sine waves

% loop over channels and trials
for chani=1:EEG.nbchan
    for triali=1:EEG.trials
        
        % hint: instantaneous frequency via interpolated random numbers
        freqmod = 20*interp1(rand(1,10),linspace(1,10,EEG.pnts));
        signal  = sin( 2*pi * ((EEG.times + cumsum(freqmod))/EEG.srate) );

        EEG.data(chani,:,triali) = signal;
    end
end

plot_simEEG(EEG,2,4)

%% 5) transient oscillations w/ Gaussian


% Gaussian and sine parameters
peaktime = 1; % seconds
width    = .12;
sinefreq = 7; % for sine wave

% create Gaussian taper
gaus = exp( -(EEG.times-peaktime).^2 / (2*width^2) );


% loop over channels and trials
for chani=1:EEG.nbchan
    for triali=1:EEG.trials
        
        % trial-unique sine wave
        cosw = cos(2*pi*sinefreq*EEG.times);
        
        EEG.data(chani,:,triali) = cosw .* gaus;
    end
end

plot_simEEG(EEG,2,5)


%% 6) repeat #3 with white noise

% list of frequencies and corresponding amplitudes
frex = [ 3 5 16 ];
amps = [ 2 4 5  ];


% loop over channels and trials
for chani=1:EEG.nbchan
    for triali=1:EEG.trials
        
        % hint: copy code from monday2
        sinewave = zeros(1,EEG.pnts);
        for si=1:numel(frex)
            sinewave = sinewave + amps(si)*sin(2*pi*frex(si)*EEG.times);
        end
        
        
        % data as a sine wave plus noise
        EEG.data(chani,:,triali) = sinewave + 5*randn(size(sinewave));
    end
end

plot_simEEG(EEG,2,6)




%% 7) repeat #5 with 1/f noise
% amount of noise
noiseamp = .3;


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
