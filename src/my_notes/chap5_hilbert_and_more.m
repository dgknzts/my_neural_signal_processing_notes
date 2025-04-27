%% VIDEO: Filter-Hilbert
%% Narrowband filtering via FIR
% filter parameters
srate   = 1024; % hz
nyquist = srate/2;
frange  = [20 25];
transw  = .1; % transition width in proportion
order   = round( 10*srate/frange(1) ); % number of time points of the filter.
% higher better but not too much. this code defines minimum useful order. 3
% is defining the three cycles.

shape   = [ 0 0 1 1 0 0 ]; % DC to nyquist. 
frex    = [ 0 frange(1)-frange(1)*transw frange frange(2)+frange(2)*transw nyquist ] / nyquist;

% filter kernel
filtkern = firls(order,frex,shape); % FIR least squared

% compute the power spectrum of the filter kernel
filtpow = abs(fft(filtkern)).^2;
% compute the frequencies vector and remove negative frequencies
hz      = linspace(0,nyquist,floor(length(filtpow)/2+ 1));
filtpow = filtpow(1:length(hz));



% plot the filter kernel
figure(20), clf
subplot(131)
plot(filtkern,'linew',2)
xlabel('Time points')
title('Filter kernel (firls)')
axis square



% plot amplitude spectrum of the filter kernel
subplot(132), hold on
plot(hz,filtpow,'ks-','linew',2,'markerfacecolor','w')
plot(frex*nyquist,shape,'ro-','linew',2,'markerfacecolor','w')


% make the plot look nicer
set(gca,'xlim',[0 frange(1)*4])
xlabel('Frequency (Hz)'), ylabel('Filter gain')
legend({'Actual';'Ideal'})
title('Frequency response of filter (firls)')
axis square


subplot(133), hold on
plot(hz,10*log10(filtpow),'ks-','linew',2,'markersize',10,'markerfacecolor','w')
plot([1 1]*frange(1),get(gca,'ylim'),'k:')
set(gca,'xlim',[0 frange(1)*4],'ylim',[-50 2])
xlabel('Frequency (Hz)'), ylabel('Filter gain (dB)')
title('Frequency response of filter (firls)')
axis square

%%% QUESTION: Is this a good filter? The answer is yes if there is a good
%             match between the "ideal" and "actual" spectral response.
%  
%%% QUESTION: One important parameter is the order (number of points in the
%             kernel). Based on your knowledge of the Fourier transform,
%             should this parameter be increased or decreased to get a
%             better filter kernel? First answer, then try it!

%% apply the filter to random noise

% generate random noise as "signal"
signal = randn(srate*4,1);

% apply the filter kernel to the signal
filtsig = filtfilt(filtkern,1,signal); % 



% <----- tangent for those without the signal-processing toolbox ----->
%  Use the following code instead of filtfilt
tmpsignal = filter(filtkern,1,signal);
tmpsignal = filter(filtkern,1,tmpsignal(end:-1:1));
filtsig1  = tmpsignal(end:-1:1);
% <----- end tangent ----->

% plot time series and its spectrum
figure(21), clf
subplot(2,4,1:3)
plot(signal,'r','linew',2)
set(gca,'xlim',[1 length(signal)])

subplot(244)
hz = linspace(0,srate,length(signal));
plot(hz,abs(fft(signal)),'r')
set(gca,'xlim',[0 frange(2)*2])


% plot time series
subplot(2,4,5:7)
plot(filtsig,'k','linew',2)
set(gca,'xlim',[1 length(signal)])
xlabel('Time (a.u.)'), ylabel('Amplitude')
title('Filtered noise in the time domain')
% plot power spectrum
subplot(248)
plot(hz,abs(fft(filtsig)) ,'k')
set(gca,'xlim',[0 frange(2)*2])
xlabel('Frequency (Hz)'), ylabel('Amplitude')
title('Spectrum of filtered signal')

%%% for comparison, plot them on top of each other
figure(22), clf
subplot(1,3,1:2), hold on
plot(signal,'r')
plot(filtsig,'k')
set(gca,'xlim',[1 length(signal)])
xlabel('Time (a.u.)'), ylabel('Amplitude')
title('Filtered noise in the time domain')
legend({'Original';'Filtered'})


subplot(133), hold on
plot(hz,abs(fft(signal)),'r')
plot(hz,abs(fft(filtsig)),'k')
set(gca,'xlim',[0 frange(2)*1.5])
xlabel('Frequency (Hz)'), ylabel('Amplitude')
title('Spectrum of filtered signal')
zoom on

%% The Hilbert transform
% take the hilbert transform
hilbFiltSig = hilbert(filtsig);

figure(23), clf
subplot(311)
plot(real(hilbFiltSig))
title('Real part of hilbert')

subplot(312)
plot(abs(hilbFiltSig))
title('Magnitude of hilbert')

subplot(313)
plot(angle(hilbFiltSig))
title('Angle of hilbert')

%% VIDEO: Short-time FFT
clear
% create signal
srate = 1000;
time  = -3:1/srate:3;
pnts  = length(time);
freqmod = exp(-time.^2)*10+10;
freqmod = freqmod + linspace(0,10,pnts);
signal  = sin( 2*pi * (time + cumsum(freqmod)/srate) );


% plot the signal
figure(24), clf
subplot(411)
plot(time,signal,'linew',1)
xlabel('Time (s)')
title('Time-domain signal')


% using MATLAB spectogram
[powspect,frex,timevec] = spectrogram(signal,hann(500),500,1000,srate);

subplot(4,1,2:4)
contourf(timevec,frex,abs(powspect),40,'linecolor','none')
set(gca,'ylim',[0 40])
colormap hot
xlabel('Time (s)'), ylabel('Frequency (Hz)')

%% a simplified version to show the mechanics

% plot the signal
figure(25), clf
subplot(411)
plot(time,signal,'linew',1)
xlabel('Time (s)')
title('Time-domain signal')

n  = 500;
nfft = 10000;
hz = linspace(0,srate,nfft);
winstep = 100;
winindx = 1:winstep:length(signal)-n;
tf = zeros(length(winindx),length(hz));
tv = zeros(length(winindx),1);

for i=1:length(winindx)
    
    % cut some signal
    datasnip = signal(winindx(i): winindx(i) + n);
    
    % compute power in this snippet
    pw = abs(fft(datasnip, nfft)).^2;


    tf(i,1:length(hz)) = pw(1:length(hz));
    
    % center time point
    tv(i) = mean(time(winindx(i): winindx(i) + n));
end

% and plot
subplot(4,1,2:4)
imagesc(tv,hz,tf')
axis xy
set(gca,'ylim',[0 40])
colormap hot
hold on 
plot(freqmod)
xlabel('Time (s)'), ylabel('Frequency (Hz)')