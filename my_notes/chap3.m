%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%% Time-domain analyses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% The theory of an ERP via simulation

% simulation details
srate   = 500; %sampling rate
time    = -1:1/srate:2; %time vector
ntrials = 100; 
sfreq   =   4.8; % sine wave frequency in Hz
gpeakt  = .43; % peak time of the frequency
gwidth  =  0.1; % gaussian width in seconds. 
% increase = lower time resolution higher freq resolution
noiseamp = 1; % noise standard deviation



% create signal
swave  = cos(2*pi*sfreq*time); % sine wave -> plot(swave)
gausw  = exp( -4*log(2)*(time-gpeakt).^2 / gwidth^2 ); % gaussian dist -> plot(gausw)
signal = swave .* gausw; %taper sine with a gaussian

% create data and multiple channels plus noise
data = repmat(signal,ntrials,1); % similar to for loop we use in chap2.
data = data + noiseamp * randn(ntrials,length(time));
% data + noise standart dev * noise

figure(1), clf
subplot(511)
plot(time,signal,'k','linew',3)
set(gca,'ylim',[-1 1])
title('Pure signal')

subplot(5,1,2:4)
imagesc(data)
set(gca,'clim',[-1 1]*noiseamp*2)
ylabel('Trial')
title('All trials: signal + noise')

subplot(515)
plot(time,mean(data),'k','linew',3)
xlabel('Time (s)')
set(gca,'ylim',[-1 1])
title('Average over trials')

%% Plotting time domain features
clear 
load('v1_laminar.mat')
% pick a channel to plot
chan2plot = 7;

figure(2), clf
subplot(311), hold on

% plot ERP from selected channel (in one line of code!)
plot(timevec, mean(csd(chan2plot,:,:),3),'b','linew',2)
set(gca,'xlim',[-.1 1.3])
title([ 'ERP from channel ' num2str(chan2plot) ])
plot([0 0],get(gca,'ylim'),'k--')
plot([0 0]+.5,get(gca,'ylim'),'k--')
plot(get(gca,'xlim'),[0 0],'k--') 
ylabel('Voltage (\muV)')


% now plot all trials from this channel
subplot(3,1,2:3)
imagesc(timevec,[],squeeze(csd(chan2plot,:,:))')
set(gca,'clim',[-1 1]*1e3,'xlim',[-.1 1.3])
xlabel('Time (s)')
ylabel('Trials')
hold on
plot([0 0],get(gca,'ylim'),'k--','linew',3)
plot([0 0]+.5,get(gca,'ylim'),'k--','linew',3)

% data for all channels

figure(3), clf

% make an image of the ERPs from all channels
contourf(timevec,1:16,squeeze(mean(csd,3)), 40, 'linecolor', 'none')
set(gca,'xlim',[-.1 1.3],'ydir','reverse')
title('Time-by-depth plot')
xlabel('Time (s)'), ylabel('Channel')
hold on
plot([0 0],get(gca,'ylim'),'k--','linew',3)
plot([0 0]+.5,get(gca,'ylim'),'k--','linew',3)

%% Filtering
% reduce data for convenience
data = double(squeeze( csd(7,:,:) ));

% cutoff frequency for low-pass filter
lowcut = 20; % in Hz

%% create and inspect the filter kernel

% filter order
filtord = round( 18 * (lowcut*1000/srate) ); 
% integer multiplied by the number of time points
% 18 is usually appropriate. higher -> better freq resolution but loosing
% time specifity. higher are resulting almost longer kernels than 
% the actual data.

% create filter
filtkern = fir1(filtord,lowcut/(srate/2),'low');

% inspect the filter kernel
% time domain
figure(1), clf
subplot(211)
plot((0:length(filtkern)-1)/srate,filtkern,'k','linew',2)
xlabel('Time (s)')
title('Time domain')

% frequency domain
subplot(212)
hz = linspace(0,srate,length(filtkern));
plot(hz,abs(fft(filtkern)).^2,'ks-','linew',2)
set(gca,'xlim',[0 lowcut*3])
xlabel('Frequency (Hz)'), ylabel('Gain')
title('Frequency domain')

%% Filtering methods
% all methods creates the almost the same output
% option 1: filter the ERP

% extract ERP
erp1 = mean(data,2);

% apply filter
erp1 = filtfilt(filtkern,1,erp1);
% filfilt() is the matlab function that applies the filtering method
% option 2: filter the single trials

erp2 = zeros(size(timevec));

for triali=1:size(data,2)
    erp2 = erp2 + filtfilt(filtkern,1,data(:,triali));
end

% complete the averaging
erp2 = erp2/triali;

% option 3: concatenate

% make one long trial
supertrial = reshape(data,1,[]);

% apply filter
supertrial = filtfilt(filtkern,1,supertrial);

% reshape back and take average
erp3 = reshape(supertrial,size(data));
erp3 = mean(erp3,2);

% plotting the before and after filters
figure(2), clf, hold on

c = 'brk';
s = 'so^';

for i=1:3
    eval([ 'plot(timevec,erp' num2str(i) ',[c(i) s(i) ''-'' ],''linew'',2)' ])
end

xlabel('Time (s)'), ylabel('Voltage (\muV)')
legend({'filter ERP';'filter trials';'filter concat'})

%% Average reference
clear
load('sampleEEGdata.mat');

% initialize new data matrices
[EEG.cardata1,EEG.cardata2] = deal( zeros(size(EEG.data)) );
% deal is helpful to initial multiple variables

% via a double-loop

for triali=1:EEG.trials
    for ti=1:EEG.pnts
        
        % new channel vector is itself minus average over channels
        EEG.cardata1(:,ti,triali) = EEG.data(:,ti,triali) - mean(EEG.data(:,ti,triali));
    end
end

% via bsxfun
% same method but different function
EEG.cardata2 = bsxfun(@minus,EEG.data,mean(EEG.data,1));
% first input: opearation
% second and third input: applying that function on those two data
% compare them 

chan2plot = 'poz';


% convert channel label to index
chanidx = strcmpi({EEG.chanlocs.labels},chan2plot);

figure(1), clf, hold on
h(1) = plot(EEG.times,mean(EEG.data(chanidx,:,:),3),'r');
h(2) = plot(EEG.times,mean(EEG.cardata1(chanidx,:,:),3),'k');
h(3) = plot(EEG.times,mean(EEG.cardata2(chanidx,:,:),3),'b');
% different references should be look the similar. this is the case here.
legend({'Earlobe';'Average (loop)';'Average (bsxfun)'})

% adjust line widths and markers
set(h,'linewidth',2)
set(h(2),'marker','s')

% adjust general plot settings
set(gca,'xlim',[-300 1200])
xlabel('Time (ms)'), ylabel('Voltage (\muV)')

%% Butterfly plot and topo-variance time series
clear
load('sampleEEGdata.mat')

figure(1), clf

% make a butterfly plot
subplot(311)
plot(EEG.times,squeeze(mean(EEG.data,3)),'linew',2)
set(gca,'xlim',[-500 1300])
title('Butterfly plot')
xlabel('Time (s)'), ylabel('Voltage (\muV)')
grid on
title('Earlobe reference')


% compute the average reference

% the fancy way...
data = bsxfun(@minus,EEG.data,mean(EEG.data,1));


subplot(312)
plot(EEG.times,squeeze(mean(data,3)),'linew',2)
set(gca,'xlim',[-500 1300])
title('Butterfly plot')
xlabel('Time (s)'), ylabel('Voltage (\muV)')
grid on
title('Average reference')


% compute the variance time series
% idea behind the plotting varience is being able to show the
% data with only single lines.

% variance for both earlobe and average reference
var_ts_ER = var( mean(EEG.data,3) );
var_ts_AR = var( mean(data,3) );

subplot(313), hold on
plot(EEG.times,var_ts_ER,'rs-','linew',2)
plot(EEG.times,var_ts_AR,'k-','linew',2)

set(gca,'xlim',[-500 1300])
xlabel('Time (s)'), ylabel('Voltage (\muV)')
grid on
title('Topographical variance time series')

%% Topography time series
clear
load('sampleEEGdata.mat')

% time points for topographies
times2plot = -200:50:800; % in ms

% convert to indices
tidx = zeros(size(times2plot));
for i=1:length(tidx)
    [~,tidx(i)] = min(abs( EEG.times-times2plot(i) ));
end

% or:
% similar to the upper method but upper is faster!
tidx = dsearchn(EEG.times',times2plot'); % they shoud be all column vectors

%topoplot time series at exact time point

figure(1), clf

% optional: redo with average reference
%EEG.data = bsxfun(@minus, EEG.data, mean(EEG.data,1));
% overall amplitudes are lower so color rescaling is neccasary


% define subplot geometry
subgeomR = ceil(sqrt(length(tidx))); % this is basically for our figure
subgeomC = ceil(length(tidx)/subgeomR); % predifining how many row and column
% we need in the figure

for i=1:length(tidx)
    subplot( subgeomR,subgeomC,i )
    topoplotIndie( mean(EEG.data(:,tidx(i),:),3),EEG.chanlocs );
    % we can also create the ERP before the loop than specify time
    % indices within that erp. which will be faster
    set(gca,'clim',[-1 1]*10)
    title([ num2str(times2plot(i)) ' ms' ])
end

% topoplot time series at average around time point
% individual time points are little bit noisy.
% so we should average in a time point window
% for example for time 50ms we choose a window between 40-60 ms and
% take average of it

% window size
twin = 10; % in ms; half of window 
% it is in ms because we take ratio of it compared to sampling rate

% convert to indices
twinidx = round(twin/(1000/EEG.srate));

figure(2), clf
for i=1:length(tidx)
    subplot( subgeomR,subgeomC,i )
    
    % time points to average together
    times2ave = tidx(i)-twinidx : tidx(i)+twinidx;
    
    % draw the topomap
    topoplotIndie( mean(mean(EEG.data(:,times2ave,:),3),2),EEG.chanlocs,'electrodes','off','numcontour',0 );
    set(gca,'clim',[-1 1]*10) % important for showing all topos in the same color limits
    title([ num2str(times2plot(i)) ' ms' ])
end


%% Important note: Some analysis methods are sensitive to some filters like
% in the chap 3 project2 solutions. peak selection differs because of the
% filter. Choosing more robust parameters are highly important 
% (e.g.,narrower time windows)