%% The discrete-time Fourier transform

% generate multi-component sine wave (exactly as before)
% define a sampling rate
srate = 1000;
% some parameters
frex   = [ 3   10   5   15   35 ];
amplit = [ 5   15   10   5   7 ];
% define time...
time = -1:1/srate:1;
% there will be some differences in the fft results compared to the ground
% truth because of the non-stationaries and number of frequencies and time
% points... will be more clear in future.

% now loop through frequencies and create sine waves
signal = zeros(1,length(time));
for fi=1:length(frex)
    signal = signal + amplit(fi)*sin(2*pi*time*frex(fi));
end

%% The Fourier transform in a loop
%   the Fourier transform in a loop, as described in lecture.
N = length(signal); % length of sequence
fourierTime = (0:N-1)/N;    % "time" used for sine waves
nyquist     = srate/2;        % Nyquist frequency -- the highest frequency you can measure in the data

% initialize Fourier output matrix
fourierCoefs = zeros(size(signal));
% These are the actual frequencies in Hz that will be returned by the
% Fourier transform. The number of unique frequencies we can measure is
% exactly 1/2 of the number of data points in the time series (plus DC).
frequencies = linspace(0,nyquist,floor(N/2)+1);

% loop over frequencies
for fi=1:N/2+1
    % create complex-valued sine wave for this frequency
    fourierSine = exp(-1i*(2*pi*fourierTime*(fi-1) + 0));
    % frequency corresp to f'th point -> f minus 1
    % compute dot product between sine wave and signal (created in the previous cell)
    fourierCoefs(fi) = dot(fourierSine,signal);
end
% scale Fourier coefficients to original scale
fourierCoefs = fourierCoefs / N;
% normalization is neccesary because dot product easily efffected by the
% number of points (longer signal = larger dot products) so we have to 
% take the average basically.

figure(7), clf
subplot(221)
plot(real(exp( -2*pi*1i*(10).*fourierTime )))
xlabel('time (a.u.)'), ylabel('Amplitude')
title('One sine wave from the FT (real part)')

subplot(222)
plot(signal)
title('Data')

subplot(212)
plot(frequencies,abs(fourierCoefs(1:length(frequencies)))*2,'*-')
xlabel('Frequency (Hz)')
ylabel('Amplitude (\muV)')
title('Amplitude spectrum derived from discrete Fourier transform')

%% Fast Fourier Transform
%%% The "slow" FT is important to understand and see implemented,
%   but in practice you should always use the FFT.
%   In this code you will see that they produce the same results.

% Compute fourier transform and scale
fourierCoefsF = fft(signal) / N;

subplot(212), hold on
plot(frequencies,abs(fourierCoefsF(1:length(frequencies)))*2,'ro')
set(gca,'xlim',[0 40])

%% Fourier coefficients are difficult to interpret 'numerically', that is,
%   it's difficult to extract the information in a Fourier coefficient
%   simply by looking at the real and imaginary parts printed out.
%   Instead, you can understand them by visualizing them (in the next cell!).
srate = 1000;
time  = (0:srate-1)/srate;
freq  = 6;

% create sine waves that differ in power and phase
sine1 = 3 * cos(2*pi*freq*time + 0 );
sine2 = 2 * cos(2*pi*freq*time + pi/2 );
sine3 = 1 * cos(2*pi*freq*time + pi*4/3 );

% compute Fourier coefficients
fCoefs1 = fft(sine1) / length(time);
fCoefs2 = fft(sine2) / length(time);
fCoefs3 = fft(sine3) / length(time);


hz = linspace(0,srate/2,floor(length(time)/2)+1);
% find the frequency of our sine wave
hz6 = dsearchn(hz',6);

% let's look at the coefficients for this frequency
disp([ '6 Hz Fourier coefficient for sin1: ' num2str(fCoefs1(hz6)) ])
disp([ '6 Hz Fourier coefficient for sin2: ' num2str(fCoefs2(hz6)) ])
disp([ '6 Hz Fourier coefficient for sin3: ' num2str(fCoefs3(hz6)) ])
% its hard to interpret in numbers so we can use polar plots.
%%% Explore the concept that the Fourier coefficients are complex numbers,
%   and can be represented on a complex plane.
%   The magnitude is the length of the line, and the phase is the angle of that line.

% make polar plots of fourier coefficients
figure(9), clf
h(1) = polarplot([0 angle(fCoefs1(hz6))],[0 2*abs(fCoefs1(hz6))],'r');
hold on
h(2) = polarplot([0 angle(fCoefs2(hz6))],[0 2*abs(fCoefs2(hz6))],'b');
h(3) = polarplot([0 angle(fCoefs3(hz6))],[0 2*abs(fCoefs3(hz6))],'m');

% adjust the plots a bit
set(h,'linewidth',5)
legend({'sine1';'sine2';'sine3'})
%% phase and power information can be extracted via Euler's formula

% extract amplitude using Pythagorian theorem
amp1 = sqrt( imag(fCoefs1).^2 + real(fCoefs1).^2 );
amp2 = sqrt( imag(fCoefs2).^2 + real(fCoefs2).^2 );
amp3 = sqrt( imag(fCoefs3).^2 + real(fCoefs3).^2 );

% extract amplitude using the Matlab function abs
% amp1 = abs( fCoefs1 );
% amp2 = abs( fCoefs2 );
% amp3 = abs( fCoefs3 );

% yet another possibility (the complex number times its conjugate)
%amp1 = fCoefs1.*conj(fCoefs1);
%%% THIS ONE IS THE FASTEST ONE

figure(9), clf

% plot amplitude spectrum
subplot(211)
plot(hz,2*amp1(1:length(hz)),'ro-','linew',3), hold on
plot(hz,2*amp2(1:length(hz)),'bp-','linew',3)
plot(hz,2*amp3(1:length(hz)),'ks-','linew',3)
set(gca,'xlim',[freq-3 freq+3])
xlabel('Frequency (Hz)'), ylabel('Amplitude')

%% and now for phase...

% extract phase angles using trigonometry
phs1 = atan(imag(fCoefs1) ./ real(fCoefs1));
phs2 = atan(imag(fCoefs2) ./ real(fCoefs3));
phs3 = atan(imag(fCoefs3) ./ real(fCoefs3));

% extract phase angles using Matlab function angle
% phs1 = angle( fCoefs1 );
% phs2 = angle( fCoefs2 );
% phs3 = angle( fCoefs3 );


% plot phase spectrum
subplot(212)
plot(hz,phs1(1:length(hz)),'ro-','linew',3), hold on
plot(hz,phs2(1:length(hz)),'bp-','linew',3)
plot(hz,phs3(1:length(hz)),'ks-','linew',3)
set(gca,'xlim',[freq-3 freq+3])
xlabel('Frequency (Hz)'), ylabel('Phase (rad.)')
% outside of the 6 hz we dont have any amplitudes so when amplitude is 0 
% phase is undefined and computer just produce random numbers for undefined
% phases. 

%% scaling of Fourier coefficients
%%% The goal of this section of code is to understand the necessity and logic
%   behind the two normalization factors that get the Fourier coefficients
%   in the same scale as the original data.


% create the signal
srate  = 1000; % hz
time   = (0:3*srate-1)/srate; % time vector in seconds
pnts   = length(time); % number of time points
signal = 2.5 * sin( 2*pi*4*time );


% prepare the Fourier transform
fourTime = (0:pnts-1)/pnts;
fCoefs   = zeros(size(signal));

for fi=1:pnts
    % create complex sine wave and compute dot product with signal
    csw = exp( -1i*2*pi*(fi-1)*fourTime );
    fCoefs(fi) = sum( signal.*csw );
end

% extract amplitudes
ampls = abs(fCoefs) / (pnts/2);
% dividing by half of the number of points is showing us the ground truth.

% compute frequencies vector
hz = linspace(0,srate/2,floor(pnts/2)+1);

figure(11), clf
stem(hz,ampls(1:length(hz)),'ks-','linew',3,'markersize',10,'markerfacecolor','w')
set(gca,'xlim',[0 10])
xlabel('Frequency (Hz)'), ylabel('Amplitude (a.u.)')
%% DC reflects the mean offset
figure(12), clf
% NOTE: below is the same signal with (1) small, (2) no, (3) large DC
%       Is the amplitude spectrum accurate at 0 Hz???
signalX1 = fft(signal+2)             / length(pnts);
signalX2 = fft(signal-mean(signal))  / length(pnts);
signalX3 = fft(signal+10)            / length(pnts);

% plot signals in the time domain
subplot(211), hold on
plot(ifft(signalX1)*pnts,'bo-')
plot(ifft(signalX2)*length(signal),'rd-')
plot(ifft(signalX3)*pnts,'k*-')
xlabel('Tme (ms)'), ylabel('Amplitude')

% plot signals in the frequency domain
subplot(212)
plot(hz,2*abs(signalX1(1:length(hz))),'bo-','linew',2,'markersize',8), hold on
plot(hz,2*abs(signalX2(1:length(hz))),'rd-','linew',2,'markersize',8)
plot(hz,2*abs(signalX3(1:length(hz))),'k*-','linew',2,'markersize',8)

xlabel('Frequencies (Hz)'), ylabel('Amplitude')
legend({'+2 mean';'de-meaned';'+10 mean'})
set(gca,'xlim',[0 10])

%% Spectral analysis of resting-state EEG
% The goal of this cell is to plot a power spectrum of resting-state EEG data.

clear
load EEGrestingState.mat
% create a time vector that starts from 0
npnts = length(eegdata);
time  = (0:npnts-1)/srate;

% plot the time-domain signal
figure(13), clf
plot(time,eegdata)
xlabel('Time (s)'), ylabel('Voltage (\muV)')
zoom on


% static spectral analysis
hz = linspace(0,srate/2,floor(npnts/2)+1);
ampl = abs(fft(eegdata)/npnts);
ampl(2:end-1) = 2*ampl(2:end-1);
powr = ampl.^2;

figure(14), clf, hold on
plot(hz,ampl(1:length(hz)),'k','linew',2)
plot(hz,powr(1:length(hz)),'r','linew',2)

xlabel('Frequency (Hz)')
ylabel('Amplitude or power')
legend({'Amplitude';'Power'})
set(gca,'xlim',[0 30])


%%% When we transform the amp into power, large values appears to be
%%% stronger and small ones weaker. Amp spectrum highlight the more subtle
%%% features, Power spectrum highlights the more prominent features.
%% Quantify alpha power over the scalp
clear
load restingstate64chans.mat
% These data comprise 63 "epochs" of resting-state data. Each epoch is a
%   2-second interval cut from ~2 minutes of resting-state.
% The goal is to compute the power spectrum of each 2-second epoch
%   separately, then average together.
%   Then, extract the average power from 8-12 Hz (the "alpha band") and
%   make a topographical map of the distribution of this power.


% convert to double-precision
EEG.data = double(EEG.data);

% two ways to get the fft for all chan and trials
% 1st way is the shortests
chanpowr = (2* abs(fft(EEG.data, [], 2) / EEG.pnts)).^2; % we didnt touch the second argument but
% we selected the third argument which is the dimension to take the fft
% over !
% 2nd way
chanpowr2 = zeros(size(EEG.data));
for chani = 1:length(EEG.data(:,1,1))
    for triali = 1:length(EEG.data(1,1,:))
        chanpowr2(chani, :, triali) = (2* abs(fft(EEG.data(chani, : , triali))/ EEG.pnts)).^2;
    end
end
% then average over trials
chanpowr = mean(chanpowr,3);

% vector of frequencies
hz = linspace(0,EEG.srate/2,floor(EEG.pnts/2)+1);
% do some plotting
% plot power spectrum of all channels
figure(15), clf
plot(hz,chanpowr(:,1:length(hz)),'linew',2)
xlabel('Frequency (Hz)'), ylabel('Power (\muV)')
set(gca,'xlim',[0 30],'ylim',[0 50])
%% now to extract alpha power
% boundaries in hz
alphabounds = [8 12];
% convert to indices
freqidx = dsearchn(hz', alphabounds');
% extract average power
alphapower = mean(chanpowr(:,freqidx(1):freqidx(2)),2);
% and plot
figure(16), clf
topoplotIndie(alphapower,EEG.chanlocs,'numcontour',0);
set(gca,'clim',[0 5])
colormap hot