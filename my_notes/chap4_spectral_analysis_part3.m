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
% SOLID FOR SEEING THE FREQUENCIES THAT ARE NOT VISIBLE..