%% Non-stationary narrowband activity via filtered noise
% we can filter (taper) the frequency domain with a gaussian to only keep a
% narrowband frequencies in the data based on the width of the gaussian.

% noise varying in time
% simulation details
pnts  = 1000;
srate =  250;

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
fc = rand(1,pnts) .* exp(1i*2*pi*rand(1,pnts)); % white noise with random phase values
% a basic random frequency spectrum

% taper with Gaussian
fc2 = fc .* fg; % we are basically multiply frequencies outside of the gaus window with 0 
% to filter them

% go back to time domain to get EEG data
signal = 2*real( ifft(fc2) );


%%% plotting
figure(1), clf
subplot(411)
plot(hz, fg)
set(gca,'xlim',[0 peakfreq*3])
xlabel('Frequency (Hz)'), ylabel('Amplitude (a.u.)')
title('Gaussian function')

subplot(412)
plot(hz, abs(fc))
set(gca,'xlim',[0 peakfreq*3])
xlabel('Frequency (Hz)'), ylabel('Amplitude (a.u.)')
title('Freq Domain Before filtering')

subplot(413)
plot(hz,abs(fc2),'k')
set(gca,'xlim',[0 peakfreq*3])
xlabel('Frequency (Hz)'), ylabel('Amplitude (a.u.)')
title('Frequency domain after filtering')

subplot(414)
plot((0:pnts-1)/srate,signal,'b')
title('Time domain')
xlabel('Time (s)'), ylabel('Amplitude')
