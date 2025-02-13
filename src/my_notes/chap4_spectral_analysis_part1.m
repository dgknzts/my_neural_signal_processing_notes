%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%% Static spectral analyses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Sine waves and parameters
% define a sampling rate
srate = 1000;
% list some frequencies
frex = [ 3   10   5   15   35 ];
% list some random amplitudes... make sure there are the same number of
% amplitudes as there are frequencies
amplit = [ 5   15   10   5   7 ];
% phases... list some random numbers between -pi and pi
phases = [  pi/7  pi/8  pi  pi/2  -pi/4 ];
% define time...
time = -1:1/srate:1; 
% odd number of time points are highly important in morelett wavelets
% to make the median time point at "0"

% now loop through frequencies and create sine waves
sine_waves = zeros(length(frex),length(time));
for fi=1:length(frex)
    sine_waves(fi,:) = amplit(fi) * sin(2*pi*time*frex(fi) + phases(fi));
end

%%% now plot
figure(1), clf
for sinei=1:length(amplit)
    subplot(length(amplit),1,sinei)
    % should be one sine wave per subplot
    plot(time,sine_waves(sinei,:),'k') %indv plots
end
figure(2), clf
subplot(211)
plot(time,sum(sine_waves,1),'k') %sum of them.
xlabel('Time (s)')

%% Check the sinewave_from_params figure. in chap4

%% Complex numbers and Euler's formula
% there are several ways to create a complex number
% be carefull about creating another variable as "i"
z = 4 + 3i;
z = 4 + 3*1i;
z = 4 + 3*sqrt(-1);
z = complex(4,3);
disp([ 'Real part is ' num2str(real(z)) ' and imaginary part is ' num2str(imag(z)) '.' ])

% plot the complex number
figure(3), clf
plot(real(z),imag(z),'s','markersize',12,'markerfacecolor','k')
% make plot look nicer
set(gca,'xlim',[-5 5],'ylim',[-5 5])
grid on, hold on, axis square
plot(get(gca,'xlim'),[0 0],'k','linew',2)
plot([0 0],get(gca,'ylim'),'k','linew',2)
xlabel('Real axis')
ylabel('Imaginary axis')
title([ 'Number (' num2str(real(z)) ' ' num2str(imag(z)) 'i) on the complex plane' ])


%% Euler's formula and the complex plane
% use Euler's formula to plot vectors
m = 4;
k = pi;
compnum = m*exp( 1i*k );

% extract magnitude and angle
mag = sqrt(compnum* conj(compnum));
phs = angle(compnum);

% plot a red dot in the complex plane corresponding to the
% real (x-axis) and imaginary (y-axis) parts of the number.
figure(4), clf
plot([0 real(compnum)],[0 imag(compnum)],'ro','linew',2,'markersize',10,'markerfacecolor','r')
% make plot look nicer
set(gca,'xlim',[-5 5],'ylim',[-5 5])
grid on, hold on, axis square
plot(get(gca,'xlim'),[0 0],'k','linew',2)
plot([0 0],get(gca,'ylim'),'k','linew',2)
xlabel('Real axis')
ylabel('Imaginary axis')
% draw a line using polar notation
h = polar([0 phs],[0 mag],'r');
set(h,'linewidth',2)
% also draw a unit circle
x = linspace(-pi,pi,100);
h = plot(cos(x),sin(x));
set(h,'color',[1 1 1]*.7) % light gray
% title
title([ 'Rectangular: [' num2str(real(compnum)) ' ' num2str(imag(compnum)) 'i ], ' ...
        'Polar: ' num2str(mag) 'e^{i' num2str(phs) '}' ])

%% dot product
%%% Create a signal (wavelet)
%   and then compute the dot product between that signal
%   and a series of sine waves.
% phase of signal
theta = 1*pi*.5;
% simulation parameters
srate = 1000;
time  = -1:1/srate:1;

% here is the signal (don't change this line)
signal = sin(2*pi*5*time + theta) .* exp( (-time.^2) / .1);
% sine wave frequencies (Hz)
sinefrex = 2:.5:10;

% plot signal
figure(4), clf
subplot(211)
plot(time,signal,'k','linew',3)
xlabel('Time (sec.)'), ylabel('Amplitude (a.u.)')
title('Signal')

dps = zeros(size(sinefrex));
amp = 1;
phase = 0;
for fi=1:length(dps)

    % create a real-valued sine wave. Note that the amplitude should be 1 and the phase should be 0
    sinew(fi,:) = amp * sin(2*pi*time*sinefrex(fi) + phase);

    % compute the dot product between sine wave and signal
    % normalize by the number of time points
    dps(:, fi) = sum(sinew(fi,:).* signal) / length(time);
end
% Very similar to fourier transform, only diff is sine waves are complex in fourier.
% and plot
subplot(212)
stem(sinefrex,dps,'k','linew',3,'markersize',10,'markerfacecolor','w')
set(gca,'xlim',[sinefrex(1)-.5 sinefrex(end)+.5], 'ylim', [-.2 .2])
xlabel('Sine wave frequency (Hz)')
ylabel('Dot product (signed magnitude)')
title([ 'Dot product of signal and sine waves (' num2str(theta) ' rad. offset)' ])

%%% Question: Try changing the 'theta' parameter. 
% What is the effect on the spectrum of dot products?
% It changes the magnitude of the signal. In different directions.
% Related to the circular angle of the pi. 

%% complex-valued sine wave

% general simulation parameters
srate = 500; % sampling rate in Hz
time  = 0:1/srate:2; % time in seconds

% sine wave parameters
freq = 5;    % frequency in Hz
ampl = 2;    % amplitude in a.u.
phas = pi; % phase in radians

% generate the complex sine wave.
csw = ampl*exp(1i*(2*pi*freq*time + phas));

% plot in 2D
figure(5), clf
subplot(211)
plot(time,real(csw), time,imag(csw),'linew',2)
xlabel('Time (sec.)'), ylabel('Amplitude')
title('Complex sine wave projections')
legend({'real';'imag'})
% plot in 3D
subplot(212)
plot3(time,real(csw),imag(csw),'k','linew',3)
xlabel('Time (sec.)'), ylabel('real part'), zlabel('imag. part')
set(gca,'ylim',[-1 1]*ampl*3,'zlim',[-1 1]*ampl*3)
axis square
rotate3d on

%% complex dot product with wavelet

% phase of signal
theta = 1*pi*3;
% simulation parameters
srate = 1000;
time  = -1:1/srate:1;

% here is the signal (don't change this line)
signal = sin(2*pi*5*time + theta) .* exp( (-time.^2) / .1);
% sine wave frequencies (Hz)
sinefrex = 2:.5:10;

% plot signal
figure(4), clf
subplot(211)
plot(time,signal,'k','linew',3)
xlabel('Time (sec.)'), ylabel('Amplitude (a.u.)')
title('Signal')

dps = zeros(size(sinefrex));
amp = 1;
phase = pi/4;
for fi=1:length(dps)

    % create a real-valued sine wave. Note that the amplitude should be 1 and the phase should be 0
    sinew(fi,:) = amp * exp(1i*2*pi*time*sinefrex(fi) + 1i*phase);

    % compute the dot product between sine wave and signal
    % normalize by the number of time points
    dps(fi) = sum(sinew(fi,:).* signal) / length(time);
end
% Very similar to fourier transform, only diff is sine waves are complex in fourier.
% and plot
subplot(212)
stem(sinefrex,abs(dps),'k','linew',3,'markersize',10,'markerfacecolor','w')
set(gca,'xlim',[sinefrex(1)-.5 sinefrex(end)+.5], 'ylim', [-.2 .2])
xlabel('Sine wave frequency (Hz)')
ylabel('Dot product (signed magnitude)')
title([ 'Dot product of signal and sine waves (' num2str(theta) ' rad. offset)' ])

%% A movie showing why complex sine waves are phase-invariant

% no need to change the code; just run and enjoy!

% create complex sine wave
csw = exp( 1i*2*pi*5*time );
rsw = cos(    2*pi*5*time );

% specify range of phase offsets for signal
phases = linspace(0,7*pi/2,100);


% setup the plot
figure(6), clf
subplot(223)
ch = plot(0,0,'ro','linew',2,'markersize',10,'markerfacecolor','r');
set(gca,'xlim',[-.2 .2],'ylim',[-.2 .2])
grid on, hold on, axis square
plot(get(gca,'xlim'),[0 0],'k','linew',2)
plot([0 0],get(gca,'ylim'),'k','linew',2)
xlabel('Cosine axis')
ylabel('Sine axis')
title('Complex plane')

% and then setup the plot for the real-dot product axis
subplot(224)
rh = plot(0,0,'ro','linew',2,'markersize',10,'markerfacecolor','r');
set(gca,'xlim',[-.2 .2],'ylim',[-.2 .2],'ytick',[])
grid on, hold on, axis square
plot(get(gca,'xlim'),[0 0],'k','linew',2)
plot([0 0],get(gca,'ylim'),'k','linew',2)
xlabel('Real axis')
title('Real number line')


for phi=1:length(phases)

    % create signal
    signal = sin(2*pi*5*time + phases(phi)) .* exp( (-time.^2) / .1);

    % compute complex dot product
    cdp = sum( signal.*csw ) / length(time);

    % compute real-valued dot product
    rdp = sum( signal.*rsw ) / length(time);

    % plot signal and real part of sine wave
    subplot(211)
    plot(time,signal, time,rsw,'linew',2)
    title('Signal and sine wave over time')

    % plot complex dot product
    subplot(223)
    set(ch,'XData',real(cdp),'YData',imag(cdp))

    % draw normal dot product
    subplot(224)
    set(rh,'XData',rdp)

    % wait a bit
    pause(.1)
end