%% Program the Fourier transform from scratch!

% random signal
N       = 50;          % length of sequence
signal  = sin(linspace(0,4*pi,N)) + randn(1,N)/2;  % data
fTime   = ((1:N)-1)/N; % "time" used in Fourier transform

% initialize Fourier output matrix
fourierCoefs = zeros(size(signal)); 

% loop over frequencies
for fi=1:N
    
    % complex sine wave for this frequency
    fourierSine = exp(-1i* 2 * pi * fTime * (fi-1));
    
    % dot product as sum of point-wise multiplications
    fourierCoefs(fi) = sum(fourierSine .* signal); %or sum(.*)
end

% divide by N to scale coefficients properly
fourierCoefs = fourierCoefs / N;

% use the fft function on the same data for comparison
fourierCoefsF = fft(signal) / N;


%% plotting

figure(1), clf

% time domain signal
subplot(211)
plot(signal)
xlabel('Time (a.u.)')
title('Data')

% amplitude of "manual" Fourier transform
subplot(212), hold on
plot(abs(fourierCoefs)*2,'*-')

% amplitude of FFT
plot(abs(fourierCoefsF)*2,'ro')


xlabel('Frequency (a.u.)')
legend({'Manual FT';'FFT'})

%% now for the IFFT

reconSig = zeros(size(signal));

for fi=1:N
    
    % create "template" sine wave
    fourierSine = fourierCoefs(fi) * exp(1i* 2 * pi * fTime * (fi-1)); %no minus sign on the inverse!
    
    % modulate by fourier coefficient and add to reconstructed signal
    reconSig = reconSig +  fourierSine; 
end
    
figure(2), clf, hold on

plot(signal,'b*-')
plot(real(reconSig),'ro')
plot(real(ifft(fourierCoefsF))*N,'ko')
xlabel('Time (a.u.)')

%% done.
%% simulation 1: frequency non-stationarity
% simulation parameters
srate = 1000;
n     = 3431;
t     = (0:n-1)/srate;

% create signals (sine wave and linear chirp)
f  = [2 10];
ff = linspace(f(1),mean(f),n);
signal1 = sin(2*pi.*ff.*t);
signal2 = sin(2*pi.*mean(f).*t);

% Fourier spectra
signal1X = fft(signal1)/n;
signal2X = fft(signal2)/n;
hz = linspace(0,srate/2,floor(n/2)+1);

% amplitude spectra
signal1Amp = abs(signal1X(1:length(hz)));
signal1Amp(2:end-1) = 2*signal1Amp(2:end-1);

signal2Amp = abs(signal2X(1:length(hz)));
signal2Amp(2:end-1) = 2*signal2Amp(2:end-1);


% plot the signals in the time domain
figure(1), clf
subplot(211)
plot(t,signal1,'b'), hold on
plot(t,signal2,'r')
xlabel('Time (sec.)'), ylabel('Amplitude')
set(gca,'ylim',[-1.1 1.1])
title('Time domain')

% and their amplitude spectra
subplot(212)
stem(hz,signal1Amp,'.-','linew',2), hold on
stem(hz,signal2Amp,'r.-','linew',2)
xlabel('Frequency (Hz)'), ylabel('Amplitude')
set(gca,'xlim',[0 20])
title('Frequency domain')
legend({'Nonstationary';'Stationary'})

%% simulation 2: amplitude non-stationarity

% simulation parameters
srate = 1000;
n     = 4000; % try other values
t     = (0:n-1)/srate;

% create signals
f = 6;
%ampmod  = 2*interp1(rand(1,10),linspace(1,10,n),'pchip');
ampmod = [ones(1,2000), zeros(1,0), ones(1,1000), zeros(1,1000)];
signal1 = ampmod .* sin(2*pi.*f.*t);
signal2 = sin(2*pi.*f.*t);

% Fourier spectra
signal1X = fft(signal1)/n;
signal2X = fft(signal2)/n;
hz = linspace(0,srate/2,floor(n/2)+1);

% amplitude spectra
signal1Amp = abs(signal1X(1:length(hz)));
signal1Amp(2:end-1) = 2*signal1Amp(2:end-1);

signal2Amp = abs(signal2X(1:length(hz)));
signal2Amp(2:end-1) = 2*signal2Amp(2:end-1);



% plot the signals in the time domain
figure(2), clf
subplot(211)
plot(t,signal1,'b'), hold on
%plot(t,signal2,'r')
xlabel('Time (sec.)'), ylabel('Amplitude')
set(gca,'ylim',[-1.1 1.1]*max(signal1))
title('Time domain')


% and their amplitude spectra
subplot(212)
stem(hz,signal1Amp,'.-','linew',2), hold on
%stem(hz,signal2Amp,'r.-','linew',2)
xlabel('Frequency (Hz)'), ylabel('Amplitude')
set(gca,'xlim',[0 20])
title('Frequency domain')
legend({'Nonstationary';'Stationary'})

%% done.

% load data
load EEGrestingState.mat
N = length(eegdata);

% time vector
timevec = (0:N-1)/srate;

% plot the data
figure(1), clf
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
figure(2), clf, hold on

plot(hz,eegpow(1:length(hz)),'k','linew',2)
plot(hzW,eegpowW/10,'r','linew',2)
set(gca,'xlim',[0 40])
xlabel('Frequency (Hz)')
legend({'"Static FFT';'Welch''s method'})

%% MATLAB pwelch (signal processing toolbox)

figure(3), clf

% create Hann window
winsize = 2*srate; % 2-second window
hannw = .5 - cos(2*pi*linspace(0,1,winsize))./2;

% number of FFT points (frequency resolution)
nfft = srate*100;

pwelch(eegdata,hannw,round(winsize/4),nfft,srate);
set(gca,'xlim',[0 40])

%% done.
