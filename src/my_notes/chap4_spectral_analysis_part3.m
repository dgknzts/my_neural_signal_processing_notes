%%Reconstruct a signal via inverse Fourier transform
% create the signal
% define a sampling rate and time vector
srate = 1000;
time  = -1:1/srate:1;
% frequencies
frex = [ 3 10 5 15 35 ];
% now loop through frequencies and create sine waves
signal = zeros(1,length(time));
for fi=1:length(frex)
    signal = signal + fi * sin(2*pi*time*frex(fi));
end

%% on to the ift!
%%% Here you will invert the Fourier transform,
%   by starting from Fourier coefficients and getting back into the time domain.

N           = length(signal); % length of sequence
fourierTime = (0:N-1)/N;    % "time" used for sine waves

reconSignal = zeros(size(signal));
fourierCoefs = fft(signal)/N;

% loop over frequencies
for fi=1:N
    % create coefficient-modulated sine wave for this frequency
    % Note: this is a complex sine wave without the minus sine in the exponential.
    fourierSine = fourierCoefs(fi) * exp(1i*(2*pi*(fi-1)*fourierTime));
    % continue building up signal...
    reconSignal = reconSignal + fourierSine;
end

% note: in practice, the inverse Fourier transform should be done using:
%reconSignal = ifft(fourierCoefs) * N;

figure(17), clf
plot(real(reconSignal),'.-')
hold on
plot(signal,'ro')
zoom on % inspect the two signals

legend({'reconstructed';'original'})

%% Frequency resolution and zero-padding
%%% We start by investigating the difference between sampling rate and
%   number of time points for Fourier frequencies.

% temporal parameters
srates  = [100 100 1000];
timedur = [  1  10    1];
freq    =     5; % in Hz
colors  = 'kmb';
symbols = 'op.';

figure(18), clf
legendText = cell(size(timedur));
for parami=1:length(colors)

    % define sampling rate in this round
    srate = srates(parami); % in Hz
    % define time
    time = -1:1/srate:timedur(parami);
    % create signal (Morlet wavelet)
    signal = cos(2*pi*freq.*time) .* exp( (-time.^2) / .05 );
    % compute FFT and normalize
    signalX = fft(signal)/length(signal);
    signalX = signalX./max(signalX);
    % define vector of frequencies in Hz
    hz = linspace(0,srate/2,floor(length(signal)/2)+1);

    % plot time-domain signal
    subplot(211)
    plot(time,signal,[colors(parami) symbols(parami) '-'],'markersize',10,'markerface',colors(parami)), hold on
    set(gca,'xlim',[-1 1])
    xlabel('Time (s)'), ylabel('Amplitude')
    title('Time domain')

    % plot frequency-domain signal
    subplot(212), hold on
    plot(hz,abs(signalX(1:length(hz))),[colors(parami) symbols(parami) '-'],'markersize',10,'markerface',colors(parami))
    xlabel('Frequency (Hz)'), ylabel('Amplitude')
    title('Frequency domain')

    legendText{parami} = [ 'srate=' num2str(srates(parami)) ', N=' num2str(timedur(parami)+1) 's' ];
end

legend(legendText)
zoom

%% zero-padding, spectral resolution, and sinc interpolation
%%% Explore the effects of zero-padding.
%   Note: I created the numbers here to show rather strong effects of zero-padding.
%    You can try with other numbers; the effects will be more mild.

figure(19), clf
signal  = [ 1 0 1 2 3 1 2 0 3 0 -1 2 1 -1];

% compute the FFT of the same signal with different DC offsets
signalX1 = fft(signal,length(signal))     / length(signal); 
% this second input is determines the length of the fft (basically padzeros automatically)
signalX2 = fft(signal, 10 + length(signal)) / length(signal); % zero-pad to 10 plus length of signal!
signalX3 = fft(signal, 100 + length(signal)) / length(signal);  % zero-pad to 100 plus length of signal!

% define frequencies vector
frex1 = linspace( 0, .5, floor(length(signalX1)/2)+1 );
frex2 = linspace( 0, .5, floor(length(signalX2)/2)+1 );
frex3 = linspace( 0, .5, floor(length(signalX3)/2)+1 );

% plot signals in the time domain
subplot(211)
plot(ifft(signalX1)*length(signal),'bo-'), hold on
plot(ifft(signalX2)*length(signal),'rd-'), hold on
plot(ifft(signalX3)*length(signal),'k*-'), hold on
xlabel('Time points (arb. units)')

% plot signals in the frequency domain
subplot(212)
plot(frex1,2*abs(signalX1(1:length(frex1))),'bo-'), hold on
plot(frex2,2*abs(signalX2(1:length(frex2))),'rd-')
plot(frex3,2*abs(signalX3(1:length(frex3))),'k*-')

xlabel('Normalized frequency units'), ylabel('Amplitude')
legend({'"Native" N';'N+10';'N+100'})
% SOLID METHOD ! FOR SEEING THE FREQUENCIES THAT ARE NOT VISIBLE..

%% Examples of sharp non-stationarities on power spectra
% sharp transitions
a = [10 2 5 8];
f = [3 1 6 12];
srate = 1000;
t = -1: 1/srate: 20;
n = length(t);
timechunks = round(linspace(1,n,length(a)+1));
% create the signal
signal = 0;
for i=1:length(a)
    signal = cat(2,signal,a(i)* sin(2*pi*f(i)*t(timechunks(i):timechunks(i+1)-1) ));
end

% compute its spectrum
signalX = fft(signal)/n;
hz = linspace(0,srate/2,floor(n/2)+1);

figure(20), clf
subplot(211)
plot(t,signal)
xlabel('Time'), ylabel('amplitude')

subplot(212)
plot(hz,2*abs(signalX(1:length(hz))),'s-','markerface','k')
xlabel('Frequency (Hz)'), ylabel('amplitude')
set(gca,'xlim',[0 20])
%% edges and edge artifacts
x = (linspace(0,1,n)>.5)+0; % +0 converts from boolean to number

% uncommenting this line shows that nonstationarities
% do not prevent stationary signals from being easily observed
x = x + .08*sin(2*pi*6*t);
% plot
figure(21), clf
subplot(211)
plot(t,x)
set(gca,'ylim',[-.1 1.1])
xlabel('Time (s.)'), ylabel('Amplitude (a.u.)')
subplot(212)
xX = fft(x)/n;
plot(hz,2*abs(xX(1:length(hz))))
set(gca,'xlim',[0 20],'ylim',[0 .1])
xlabel('Frequency (Hz)'), ylabel('Amplitude (a.u.)')
%% Examples of smooth non-stationarities on power spectra
srate = 1000;
t = 0:1/srate:10;
n = length(t);
f = 3; % frequency in Hz

% sine wave with time-increasing amplitude
ampl1 = linspace(1,10,n);
%ampl1 = abs(interp1(linspace(t(1),t(end),10),10*rand(1,10),t,'spline'));
ampl2 = mean(ampl1);

signal1 = ampl1 .* sin(2*pi*f*t);
signal2 = ampl2 .* sin(2*pi*f*t);

% obtain Fourier coefficients and Hz vector
signal1X = fft( signal1 )/n;
signal2X = fft( signal2 )/n;
hz = linspace(0,srate,n);

figure(22), clf
subplot(211)
plot(t,signal2,'r','linew',2), hold on
plot(t,signal1,'linew',2)
xlabel('Time'), ylabel('amplitude')
subplot(212)
plot(hz,2*abs(signal2X(1:length(hz))),'ro-','markerface','r','linew',2)
hold on
plot(hz,2*abs(signal1X),'s-','markerface','k','linew',2)
xlabel('Frequency (Hz)'), ylabel('amplitude')
set(gca,'xlim',[1 7])
legend({'Stationary';'Non-stationary'})

%% frequency non-stationarity
f  = [2 10];
ff = linspace(f(1),mean(f),n);
signal1 = sin(2*pi.*ff.*t);
signal2 = sin(2*pi.*mean(ff).*t);

signal1X = fft(signal1)/n;
signal2X = fft(signal2)/n;
hz = linspace(0,srate/2,floor(n/2));

figure(23), clf

subplot(211)
plot(t,signal1), hold on
plot(t,signal2,'r')
xlabel('Time'), ylabel('amplitude')
set(gca,'ylim',[-1.1 1.1])

subplot(212)
plot(hz,2*abs(signal1X(1:length(hz))),'.-'), hold on
plot(hz,2*abs(signal2X(1:length(hz))),'r.-')
xlabel('Frequency (Hz)'), ylabel('amplitude')
set(gca,'xlim',[0 20])
%% examples of rhythmic non-sinusoidal time series
% parameters
srate = 1000;
time  = (0:srate*6-1)/srate;
npnts = length(time);
hz = linspace(0,srate,npnts);

% various mildly interesting signals to test
signal = detrend( sin( cos(2*pi*time)-1 ) );
%signal = sin( cos(2*pi*time) + time );
%signal = detrend( cos( sin(2*pi*time).^4 ) );


% plot!
figure(24), clf

% time domain
subplot(211)
plot(time,signal,'k','linew',3)
xlabel('Time (s)')

% frequency domain
subplot(212)
plot(hz,abs(fft(signal)),'k','linew',3)
set(gca,'xlim',[0 20])
ylabel('Frequency (Hz)')

%% Welch's method on phase-slip data
srate = 1000;
time  = (0:srate-1)/srate;

signal = [ sin(2*pi*10*time) sin(2*pi*10*time(end:-1:1)) ];

figure(25), clf
subplot(211)
plot(signal)

subplot(223)
bar(linspace(0,srate,length(signal)),2*abs(fft(signal))/length(signal))
set(gca,'xlim',[5 15])
title('Static FFT')
xlabel('Frequency (Hz)')

% Now for Welch's method
% parameters
winlen = 500; % window length in points (same as ms if srate=1000!)
skip = 100; % also in time points;
% vector of frequencies for the small windows
hzL = linspace(0,srate/2,floor(winlen/2)+1);
nwindow = length(time) / winlen;
timesize = (nwindow + (nwindow -1));
% initialize time-frequency matrix
welchspect = zeros(1, length(hzL));
% Hann taper
hwin = .5*(1-cos(2*pi*(1:winlen) / (winlen-1)));

% loop over time windows
nbins = 0; % note: in the video this is nbins=1 but that's a typo.
for ti=1:skip:length(signal)-winlen
    
    % extract part of the signal
    tidx    = ti:ti+winlen-1;
    tmpdata = signal(tidx);
    
    % FFT of these data (does the taper help?)
    x = fft(hwin.*tmpdata) / winlen;
    
    % and put in matrix
    welchspect = welchspect + 2*abs(x(1:length(hzL)));
    nbins = nbins + 1;
end

% divide by nbins to complete average
welchspect = welchspect/nbins;

subplot(224)
bar(hzL,welchspect)
set(gca,'xlim',[5 15])
title('Welch''s method')
xlabel('Frequency (Hz)')
% Freq res highly reduced..



%% Welch's method on resting-state EEG data
% load data
addpath 'G:\My Drive\signal_processing_mike_cohen\data\chap4'
load EEGrestingState.mat
N = length(eegdata);

% time vector
timevec = (0:N-1)/srate;

% plot the data
figure(26), clf
plot(timevec,eegdata,'k')
xlabel('Time (seconds)'), ylabel('Voltage (\muV)')

%% one big FFT (not Welch's method)

% "static" FFT over entire period, for comparison with Welch
eegpow = abs( fft(eegdata)/N ).^2;
hz = linspace(0,srate/2,floor(N/2)+1);

%% "manual" Welch's method

% window length in seconds*srate
winlength = 1*srate;
% number of points of overlap
nOverlap = round(srate/2);

% NOTE about the variable nOverlap: 
% This variable actually defines the number of data 
% points to skip forwards. The true number of overlapping 
% points is winlength-nSkip. Apologies for the confusion,
% and thanks to Eleonora De Filippi for catching that mistake.

% window onset times
winonsets = 1:nOverlap:N-winlength;

% note: different-length signal needs a different-length Hz vector
hzW = linspace(0,srate/2,floor(winlength/2)+1);

% Hann window
hannw = .5 - cos(2*pi*linspace(0,1,winlength))./2;

% initialize the power matrix (windows x frequencies)
eegpowW = zeros(1,length(hzW));

% loop over frequencies
for wi=1:length(winonsets)
    
    % get a chunk of data from this time window
    datachunk = eegdata(winonsets(wi):winonsets(wi)+winlength-1);
    
    % apply Hann taper to data
    datachunk = datachunk .* hannw;
    
    % compute its power
    tmppow = abs(fft(datachunk)/winlength).^2;
    
    % enter into matrix
    eegpowW = eegpowW  + tmppow(1:length(hzW));
end

% divide by N
eegpowW = eegpowW / length(winonsets);


%%% plotting
figure(27), clf, subplot(211), hold on

plot(hz,eegpow(1:length(hz)),'k','linew',2)
plot(hzW,eegpowW/10,'r','linew',2)
set(gca,'xlim',[0 40])
xlabel('Frequency (Hz)')
legend({'"Static FFT';'Welch''s method'})
title('Using FFT and Welch''s')

%% MATLAB pwelch

subplot(212)

% create Hann window
winsize = 2*srate; % 2-second window
hannw = .5 - cos(2*pi*linspace(0,1,winsize))./2;

% number of FFT points (frequency resolution)
nfft = srate*100;

pwelch(eegdata,hannw,round(winsize/4),nfft,srate);
set(gca,'xlim',[0 40])


%% VIDEO: Welch's method on v1 laminar data
clear
load v1_laminar.mat
csd = double(csd);

% specify a channel for the analyses
chan2use = 7;
% create Hann window
hannw = .5 - cos(2*pi*linspace(0,1,size(csd,2)))./2;
% Welch's method using MATLAB pwelch
[pxx,hz] = pwelch(squeeze(csd(chan2use,:,:)),hannw,round(size(csd,2)/10),1000,srate);

figure(28), clf
subplot(211)
plot(timevec,mean(csd(chan2use,:,:),3),'linew',2)
set(gca,'xlim',timevec([1 end]))
xlabel('Time (s)'), ylabel('Voltage (\muV)')

subplot(212)
plot(hz,mean(pxx,2),'linew',2)
set(gca,'xlim',[0 140])
xlabel('Frequency (Hz)')
ylabel('Power (\muV^2)')

%% done.