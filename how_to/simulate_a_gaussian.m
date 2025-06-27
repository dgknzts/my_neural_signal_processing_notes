%% Gaussian
% formula 1: amplitude * e^ (- (time - peak time)^2/ 2 s^2)
% amplitude : strength ? y axis.
% e : exponential : exp() function
% peak time : time point at the peak
% s : full width at half maximum.
% formula 2 : amplitude * exp ((-4 * log(2) * time^2) / h ^2)
% h = full width at half maximum (fwhm) in seconds.
% second formula is easier for expressing the main parameter "h".
% log part works as a normalization factor.

% simulation parameters
time  = -3:1/1000:3;
ptime = -1.3131;  % peak time 
ampl  = 30; % amplitude
fwhm  = 1;

% Gaussian
gwin = ampl* exp(-(4*log(2)*(time-ptime).^2) / fwhm^2);

% empirical FWHM
gwinN   = gwin./max(gwin); %normalizing the amp to 0 to 1. nec when determining the fwhm
midp    = dsearchn(time',ptime); % finding the peak time point
pst5    = midp-1+dsearchn(gwinN(midp:end)',.5); % pre half maximum time point
pre5    = dsearchn(gwinN(1:midp)',.5); % post half maximum time point
empfwhm = time(pst5) - time(pre5); %full width at half maximum


% visualize
figure(3), clf, hold on
plot(time,gwin,'k','linew',2)
plot(time([pre5 pst5]),gwin([pre5 pst5]),'ro--','markerfacecolor','k')
plot(time([pre5 pre5]),[0 gwin(pre5)],'r:')
plot(time([pst5 pst5]),[0 gwin(pst5)],'r:')
title([ 'Requested FWHM: ' num2str(fwhm) 's, empirical FWHM: ' num2str(empfwhm) 's' ])
xlabel('Time (s)'), ylabel('Amplitude')