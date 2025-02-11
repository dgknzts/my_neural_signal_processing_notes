%  VIDEO: ISPC (Inter-site phase clustering)

addpath("G:\My Drive\signal_processing_mike_cohen\src\data\synchronization")
% load in V1 mouse data
load v1_laminar

%% setup connectivity parameters

% The goal here is to explore the mechanism of phase synchronization.
% extract phase angle time series from two channels from trial 1 at 8 Hz.


% channels for connectivity
chan1idx = 1;
chan2idx = 8;


% create a complex Morlet wavelet (don't forget to plot the wavelet!!!)
cent_freq = 8;
time      = -1.5:1/srate:1.5;
s         = 8/(2*pi*cent_freq);
wavelet   = exp(2*1i*pi*cent_freq.*time) .* exp(-time.^2./(2*s^2));
half_wavN = (length(time)-1)/2;

% FFT parameters
nWave = length(time);
nData = size(csd,2);
nConv = nWave + nData - 1;

% FFT of wavelet (check nfft)
waveletX = fft(wavelet,nConv);
waveletX = waveletX ./ max(waveletX);

% initialize output time-frequency data
phase_data = zeros(2, size(csd,2));
real_data  = zeros(2, size(csd,2));



% analytic signal of channel 1
dataX = fft(mean(csd(chan1idx,:,:),3), nConv, 2);
as = ifft(waveletX.*dataX,nConv);
as = as(half_wavN+1:end-half_wavN);

% collect real and phase data
phase_data(1,:) = angle(as); % extract phase angles
real_data(1,:)  = real(as);  % extract the real part (projection onto real axis)

% analytic signal of channel 2
dataX = fft(mean(csd(chan2idx,:,:),3), nConv, 2);
as = ifft(waveletX.*dataX,nConv);
as = as(half_wavN+1:end-half_wavN);

% collect real and phase data
phase_data(2,:) = angle(as);
real_data(2,:)  = real(as);

%% setup figure and define plot handles

% note: This cell is just setting up the figure for the following cell. 
%       You can run it and move on.


% open and name figure
figure(1), clf
set(gcf,'NumberTitle','off','Name','Movie magic minimizes the magic.');

% draw the filtered signals
subplot(221)
filterplotH1 = plot(timevec(1,:),real_data(1,:),'b');
hold on
filterplotH2 = plot(timevec(1,:),real_data(2,:),'m');
set(gca,'xlim',[timevec(1) timevec(end)],'ylim',[min(real_data(:)) max(real_data(:))])
xlabel('Time (ms)')
ylabel('Voltage (\muV)')
title([ 'Filtered signal at ' num2str(cent_freq) ' Hz' ])

% draw the phase angle time series
subplot(222)
phaseanglesH1 = plot(timevec(1,:),phase_data(1,:),'b');
hold on
phaseanglesH2 = plot(timevec(1,:),phase_data(2,:),'m');
set(gca,'xlim',[timevec(1) timevec(end)],'ylim',[-pi pi]*1.1)
xlabel('Time (ms)')
ylabel('Phase angle (radian)')
title('Phase angle time series')

% draw phase angles in polar space
subplot(223)
polar2chanH1 = polar([zeros(1,1) phase_data(1,1)]',repmat([0 1],1,1)','b');
hold on
polar2chanH2 = polar([zeros(1,1) phase_data(2,1)]',repmat([0 1],1,1)','m');
title('Phase angles from two channels')

% draw phase angle differences in polar space
subplot(224)
polarAngleDiffH = polar([zeros(1,1) phase_data(2,1)-phase_data(1,1)]',repmat([0 1],1,1)','k');
title('Phase angle differences from two channels')

%% now update plots at each timestep

for ti=1:5:length(timevec)
    
    % update filtered signals
    set(filterplotH1,'XData',timevec(1:ti),'YData',real_data(1,1:ti))
    set(filterplotH2,'XData',timevec(1:ti),'YData',real_data(2,1:ti))
    
    % update cartesian plot of phase angles
    set(phaseanglesH1,'XData',timevec(1:ti),'YData',phase_data(1,1:ti))
    set(phaseanglesH2,'XData',timevec(1:ti),'YData',phase_data(2,1:ti))
    
    subplot(223), cla
    polar([zeros(1,ti) phase_data(1,1:ti)]',repmat([0 1],1,ti)','b');
    hold on
    polar([zeros(1,ti) phase_data(2,1:ti)]',repmat([0 1],1,ti)','m');
    
    subplot(224), cla
    polar([zeros(1,ti) phase_data(2,1:ti)-phase_data(1,1:ti)]',repmat([0 1],1,ti)','k');
    
    drawnow
end

%% now quantify phase synchronization between the two channels

% phase angle differences
phase_angle_differences = phase_data(2, : ) - phase_data(1,:);

% euler representation of angles
euler_phase_differences = exp(1i*phase_angle_differences);

% mean vector (in complex space)
mean_complex_vector = mean(euler_phase_differences);

% length of mean vector (this is the "M" from Me^ik, and is the measure of phase synchronization)
phase_synchronization = abs(mean_complex_vector);

disp([ 'Synchronization between ' num2str(chan1idx) ' and ' num2str(chan2idx) ' is ' num2str(phase_synchronization) '!' ])

% of course, this could all be done on one line:
phase_synchronization = abs(mean(exp(1i*phase_data(2, : ) - phase_data(1,:))));

% notice that the order of subtraction is meaningless (see below), which means that this measure of synchronization is non-directional!
phase_synchronization_backwards = abs(mean(exp(1i*phase_data(1, : ) - phase_data(2,:))));


% now plot mean vector
subplot(224)
hold on
h = polar([0 angle(mean_complex_vector)],[0 phase_synchronization]);
set(h,'linewidth',6,'color','g')


%%  VIDEO: Laplacian in simulated EEG data
% The goal here is to build intuition for the Laplacian as a high-pass
% spatial filter using simulated EEG data.

clear
load sampleEEGdata.mat

chan1 = 'pz';
chan2 = 'cz';

chan1idx = strcmpi(chan1,{EEG.chanlocs.labels});
chan2idx = strcmpi(chan2,{EEG.chanlocs.labels});


% compute inter-channel distances
eucdist = zeros(2,64);
for chani=1:EEG.nbchan
    eucdist(1,chani) = sqrt( (EEG.chanlocs(chani).X-EEG.chanlocs(chan1idx).X)^2 + (EEG.chanlocs(chani).Y-EEG.chanlocs(chan1idx).Y)^2 + (EEG.chanlocs(chani).Z-EEG.chanlocs(chan1idx).Z)^2 );
    eucdist(2,chani) = sqrt( (EEG.chanlocs(chani).X-EEG.chanlocs(chan2idx).X)^2 + (EEG.chanlocs(chani).Y-EEG.chanlocs(chan2idx).Y)^2 + (EEG.chanlocs(chani).Z-EEG.chanlocs(chan2idx).Z)^2 );
end

% generate low- and high- spatial frequency activity
lo_spatfreq = 50*exp(- (eucdist(1,:).^2)/ 300000 ); 
hi_spatfreq =    exp(- (eucdist(2,:).^2)/   3000 );


% compute Laplacian of summed activity
surf_lap = laplacian_perrinX(hi_spatfreq+lo_spatfreq,[EEG.chanlocs.X],[EEG.chanlocs.Y],[EEG.chanlocs.Z]);


% and plot
figure(2), clf
subplot(221)
topoplotIndie(lo_spatfreq,EEG.chanlocs,'numcontour',0);
title('Low spatial frequency feature')

subplot(222)
topoplotIndie(hi_spatfreq,EEG.chanlocs,'numcontour',0);
title('High spatial frequency features')

subplot(223)
topoplotIndie((lo_spatfreq + hi_spatfreq),EEG.chanlocs,'numcontour',0);
title('Low+high features')

subplot(224)
topoplotIndie(surf_lap,EEG.chanlocs,'numcontour',0);
title('Laplacian of low+high features')



%% VIDEO: Laplacian in real EEG data
% Now we will apply the Laplacian to real EEG data. We'll look at the
% topographical maps and ERPs before vs. after the Laplacian.
clear
load sampleEEGdata.mat

% pick a time point for the topoplot, and pick a channel for the ERPs
time2plot = 250; % ms
chan2plot = 'p3';


% convert to indices
tidx = dsearchn(EEG.times',time2plot);
chanidx = strcmpi({EEG.chanlocs.labels},chan2plot);


% Compute the laplacian and store as a new field in the EEG structure.
EEG.lap = laplacian_perrinX(squeeze(EEG.data(:, :, :)), [EEG.chanlocs.X],[EEG.chanlocs.Y],[EEG.chanlocs.Z]);


% The voltage and Laplacian data are in different scales. To compare them
% directly, they need to be independently normalized (z-transform).
voltERP = mean(EEG.data(chanidx,:,:),3);
voltERP = (voltERP - mean(voltERP)) / std(voltERP);

lapERP = mean(EEG.lap(chanidx,:,:),3);
lapERP = (lapERP - mean(lapERP)) / std(lapERP);



% now plot
figure(5), clf
subplot(221)
topoplotIndie(mean(EEG.data(:,tidx,:),3),EEG.chanlocs,'electrodes','labels','numcontour',0);
title([ 'Voltage (' num2str(time2plot) ' ms)' ])

subplot(222)
topoplotIndie(mean(EEG.lap(:,tidx,:),3),EEG.chanlocs,'electrodes','numbers','numcontour',0);
title([ 'Laplacian (' num2str(time2plot) ' ms)' ])

subplot(212)
plot(EEG.times,voltERP, EEG.times,lapERP,'linew',2)
set(gca,'xlim',[-300 1200])
legend({'Voltage';'Laplacian'})
title([ 'ERP from channel ' chan2plot ])
xlabel('Time (ms)'), ylabel('Data (z-score)')

%%% QUESTION:
%   Based on looking at the topoplots, pick an electrode that you think
%   will look (1) very different, and (2) very similar, for voltage and
%   Laplacian. Then show ERPs from those channels. Why do time courses 
%   from some channels look more similar than time courses from other
%   channels?