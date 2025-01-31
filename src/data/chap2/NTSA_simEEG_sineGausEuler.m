%%
%     COURSE: Solved challenges in neural time series analysis
%    SECTION: Simulating EEG data
%      VIDEO: The three important equations (sine, Gaussian, Euler's)
% Instructor: sincxpress.com
%
%%


%% intuition about sine waves

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

%% the sum of sine waves can appear like a complicated time series

% create a few sine waves and sum

% define a sampling rate
srate = 1000;

% list some frequencies
frex = [ 4.8   6   7.5 ];

% list some random amplitudes... make sure there are the same number of
% amplitudes as there are frequencies!
amplit = [5 2 5 ];

% phases... list some random numbers between -pi and pi
phases = [  pi/7  pi/8  pi  pi/2  -pi/4 ];

% define time...
time = -1:1/srate:1;


% now we loop through frequencies and create sine waves
sine_waves = zeros(length(frex),length(time));
for fi=1:length(frex)
    sine_waves(fi,:) = amplit(fi) * sin(2*pi*time*frex(fi) + phases(fi));
end

% now plot the result
figure(2), clf
plot(time,sum(sine_waves))
title('sum of sine waves')
xlabel('Time (s)'), ylabel('Amplitude (arb. units)')


% now plot each wave separately
figure(3), clf
for fi=1:length(frex)
    subplot(length(frex),1,fi)
    plot(time, sine_waves(fi,:))
    axis([ time([1 end]) -max(amplit) max(amplit) ])
end

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

%% done.
