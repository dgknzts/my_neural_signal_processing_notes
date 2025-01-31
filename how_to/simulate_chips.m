%% chirps: simulated sine waves varying frequency acrosss time
% signal will transition linearly from first frequency to the second
%   frequency

% formula: sin(2*pi* (y(t) + t(t)) -> this is for only generating one time point at time t.
%   y and t are vectors that are functions of time. they are changing in
%   every time point.

% y(t) = delta  *  cumsum( vector of frequencies ) across all time points.

% delta : sampling interval -> inverse of the sampling rate. e.g.: if sampling
%   rate is 1000 hz delta will be 1. distance between time points. if
%   sampling rate is 250 hz delta will be 4. 

% simulation details
pnts  = 2000;
srate = 1024;
time  = (0:pnts-1)/srate;

% "bipolar" chirp 
freqmod = linspace(5,15,pnts);

% multipolar chirp / uncomment this to see the difference between a linear
% frequency range
%k = 10; % poles for frequencies
%freqmod = 20*interp1(rand(1,k),linspace(1,k,pnts)); % creates k random
%numbers and linearly interpolate them. 

% signal time series
signal  = sin( 2*pi * ((time + cumsum(freqmod))/srate) );

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