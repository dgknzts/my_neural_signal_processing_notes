%% Getting to know Morlet wavelets
% parameters
srate = 1000;         % in hz
time  = -1:1/srate:1; % best practice is to have time=0 at the center of the wavelet
frex  = 4.8;         % frequency of wavelet, in Hz

% create sine wave (actually cosine, just to make it nice and symmetric)
sine_wave = cos(2 * pi * time * frex);

% create Gaussian window
fwhm = .5; % width of the Gaussian in seconds
gaus_win = exp( (-4*log(2)*time.^2) / (fwhm^2) );

% now create Morlet wavelet
mw = sine_wave .* gaus_win;

figure(1), clf
subplot(211), hold on
plot(time,sine_wave,'r')
plot(time,gaus_win,'b')
plot(time,mw,'k','linew',3)
xlabel('Time (s)'), ylabel('Amplitude')
legend({'Sine wave';'Gaussian';'Morlet wavelet'})
title('Morlet wavelet in the time domain')


% Morlet wavelet in the frequency domain
% confirm that the shape of the power spectrum of a Morlet wavelet is gaussian
pnts = length(time);
mwX = 2*abs(fft( mw )/pnts); % uh oh
hz  = linspace(0,srate,pnts);

subplot(212)
plot(hz,mwX,'k','linew',2)
xlabel('Frequency (Hz)')
ylabel('Amplitude')
title('Morlet wavelet in the frequency domain')

% Observations: - Notice that the amplitude spectrum is symmetric.
%               - Also notice the peak amplitude in the 
%                 time vs. frequency domains.
%               - The Hz units are incorrect above Nyquist. This is just a
%                 convenience plotting trick.
% 
%   TO DO: Change the following parameters to observe the effects:
%          - frex
%          - fwhm
%          - time (start and end values)

%% Time-domain convolution
% first example to build intuition

% make a kernel (e.g., Gaussian)
kernel = exp( -linspace(-2,2,20).^2 );
kernel = kernel./sum(kernel);

% try these options
% kernel = -kernel;
% kernel = [ zeros(1,9) 1 -1 zeros(1,9) ]; % edge detector!

% create the signal
signal = [ zeros(1,30) ones(1,2) zeros(1,20) ones(1,30) 2*ones(1,10) zeros(1,30) -ones(1,10) zeros(1,40) ];

% plot!
figure(2), clf, hold on
plot(kernel+1,'b','linew',3)
plot(signal,'k','linew',3)
% use MATLAB conv function for now
plot( conv(signal,kernel, 'same') ,'r','linew',3)

set(gca,'xlim',[0 length(signal)])
legend({'thing';'stuff';'Roccinante'})

%% a simpler example in more detail
% create signal
signal = zeros(1,20);
signal(8:15) = 1;

% create convolution kernel
kernel = [.6 .4 .2 0 -0.2];
%kernel = kernel - mean(kernel);
% convolution sizes
nSign = length(signal);
nKern = length(kernel);
nConv = nSign + nKern - 1;


figure(3), clf
% plot the signal
subplot(311)
plot(signal,'o-','linew',2,'markerface','g','markersize',9)
%set(gca,'ylim',[-.1 1.1],'xlim',[1 nSign])
title('Signal')

% plot the kernel
subplot(312)
plot(kernel,'o-','linew',2,'markerface','r','markersize',9)
%set(gca,'xlim',[1 nSign],'ylim',[-.5 .5])
title('Kernel')


% plot the result of convolution
subplot(313)
plot(conv(signal,kernel,'same'),'o-','linew',2,'markerface','b','markersize',9)
%set(gca,'xlim',[1 nSign],'ylim',[-.6 1])
title('Result of convolution')

%% convolution in animation
%%% Just run the whole cell and enjoy!
% movie time parameter!
refresh_speed = .6; % seconds

half_kern = floor(nKern/2);

% flipped version of kernel
kflip = kernel(end:-1:1);%-mean(kernel);

% zero-padded data for convolution
dat4conv = [ zeros(1,half_kern) signal zeros(1,half_kern) ];

% initialize convolution output
conv_res = zeros(1,nConv);

%%% initialize plot
figure(4), clf, hold on
plot(dat4conv,'o-','linew',2,'markerface','g','markersize',9)
hkern = plot(kernel,'o-','linew',2,'markerface','r','markersize',9);
hcres = plot(kernel,'s-','linew',2,'markerface','k','markersize',15);
set(gca,'ylim',[-1 1]*3,'xlim',[0 nConv+1])
plot([1 1]*(half_kern+1),get(gca,'ylim'),'k--')
plot([1 1]*(nConv-2),get(gca,'ylim'),'k--')
legend({'Signal';'Kernel (flip)';'Convolution'})

% run convolution
for ti=half_kern+1:nConv-half_kern
    
    % get a chunk of data
    tempdata = dat4conv(ti-half_kern:ti+half_kern);
    
    % compute dot product (don't forget to flip the kernel backwards!)
    conv_res(ti) = sum( tempdata.*kflip );
    
    % update plot
    set(hkern,'XData',ti-half_kern:ti+half_kern,'YData',kflip);
    set(hcres,'XData',half_kern+1:ti,'YData',conv_res(half_kern+1:ti))
    
    pause(refresh_speed)
end
%% The five steps of convolution
% now for convolution via spectral multiplication

% Step 1: N's of convolution
ndata = length(signal);
nkern = length(kernel);
nConv = ndata+nkern - 1;% length of result of convolution
halfK = floor(nkern/2);

% Step 2: FFTs
dataX = fft( signal, nConv); % important: make sure to properly zero-pad!
kernX = fft( kernel, nConv);

% Step 3: multiply spectra
convresX = dataX .* kernX;

% Step 4: IFFT
convres = ifft(convresX);

% Step 5: cut off "wings"
convres = convres(halfK+1:end-halfK+1);

%%% and plot for confirmation!
plot( convres,'go','markerfacecolor','g' )


%%  VIDEO: Convolve real data with a Gaussian
clear

load v1_laminar.mat

% signal will be ERP from channel 7
signal = mean(csd(7,:,:),3);

% create a Gaussian
h = 0.01; % FWHM in seconds

gtime = -1:1/srate:1;
gaus = exp( -4*log(2)*gtime.^2 / h^2 );
gaus = gaus./sum(gaus); % amplitude normalization


%%%% run convolution
% Step 1: N's of convolution
ndata = length(signal);
nkern = length(gaus);
nConv = ndata+nkern - 1;% length of result of convolution
halfK = floor(nkern/2);

% Step 2: FFTs
dataX = fft( signal, nConv ); % important: make sure to properly zero-pad!
kernX = fft( gaus, nConv );

% Step 3: multiply spectra
mult = dataX .* kernX;

% Step 4: inverse FFT to get back to the time domain
convres  = ifft( mult );

% Step 5: cut off "wings"
convres = convres(halfK+1:end-halfK+1);


% plotting!
figure(6), clf, hold on
plot(timevec,signal)
plot(timevec,convres,'r','linew',2)

set(gca,'xlim',[-.1 1.4])
legend({'Original ERP';'Gaussian-convolved'})
xlabel('Time (s)'), ylabel('Activity (\muV)')

% show the mechanism of convolution (spectral multiplication)

hz = linspace(0,srate,nConv);

figure(7), clf, hold on
plot(hz,abs(dataX)./max(abs(dataX))); % normalized for visualization
plot(hz,abs(kernX));
plot(hz,abs(dataX.*kernX)./max(abs(dataX)),'k','linew',2)

set(gca,'xlim',[0 150])
xlabel('Frequency (Hz)'), ylabel('Amplitude (norm.)')
legend({'Original signal';'Kernel';'Convolution result'})

%%  VIDEO: Complex Morlet wavelets
% setup parameters
srate = 1000;         % in hz
time  = -1:1/srate:1; % best practice is to have time=0 at the center of the wavelet
frex  = 2*pi;         % frequency of wavelet, in Hz

% create sine wave
sine_wave = exp( 1i * 2 *pi * frex * time ); % hmmmm

% create Gaussian window
fwhm = .5; % width of the Gaussian in seconds
gaus_win = exp( (-4* log(2) * time.^2 )/ fwhm^2 );

% now create Morlet wavelet
cmw = sine_wave .* gaus_win;

figure(8), clf

subplot(211), hold on
plot(time_wavelet,real(cmw),'b')
plot(time_wavelet,imag(cmw),'r--')
xlabel('Time (s)'), ylabel('Amplitude')
legend({'real part';'imag part'})
title('Complex Morlet wavelet in the time domain')

% complex Morlet wavelet in the frequency domain

pnts = length(time);

mwX = abs(fft( cmw)/EEG.pnts);
hz  = linspace(0,EEG.srate,EEG.pnts);

subplot(212)
plot(hz,mwX,'k','linew',2)
xlabel('Frequency (Hz)')
ylabel('Amplitude')
title('Complex Morlet wavelet in the frequency domain')
%% Complex Morlet wavelet convolution
clear

load v1_laminar.mat

% extract a bit of data for convenience
data = csd(6,:,10) ;


% create a complex Morlet wavelet
time = (0:2*srate)/srate;
time = time - mean(time); % note the alternative method for creating centered time vector 
frex = 45; % frequency of wavelet, in Hz

% create Gaussian window
s = 7 / (2*pi*frex); % using num-cycles formula
cmw  = exp(1i*2*pi*frex*time) .* exp( -time.^2/(2*s^2) );


%%% now for convolution

% Step 1: N's of convolution
nSignal = length(data);
nKern = length(cmw);
nConv = nSignal + nKern - 1;
halfK = floor(nKern / 2);

% Step 2: FFTs
kernX = fft(cmw, nConv);
signalX = fft(data, nConv);

% Step 2.5: normalize the wavelet (try it without this step!)
kernX = kernX ./ max(kernX);

% Step 3: multiply spectra
convresX = signalX .* kernX;

% Step 4: IFFT
convres = ifft(convresX);

% Step 5: cut off "wings"
convres = convres(halfK+1:end-halfK+1);

% plotting
% compute hz for plotting
hz = linspace(0,srate/2,floor(nConv/2)+1);

figure(9), clf
subplot(211),  hold on

% plot power spectrum of data
plot(hz,abs(signalX(1:length(hz))),'b')

% plot power spectrum of wavelet
plot(hz,abs(kernX(1:length(hz))).*max(abs(signalX))/2)

% plot power spectrum of convolution result
plot(hz,abs(convresX(1:length(hz))),'k','linew',2)
set(gca,'xlim',[0 frex*2])
xlabel('Frequency (Hz)'), ylabel('Amplitude (a.u.)')
legend({'Data spectrum';'Wavelet spectrum';'Convolution result'})


%%% now plot in the time domain
subplot(212), hold on
plot(timevec,data,'b')
plot(timevec,real(convres),'k','linew',2)
set(gca,'xlim',[-.1 1.3])
legend({'LFP data';'Convolution result'})
xlabel('Time (s)'), ylabel('Activity (\muV)')

%% extracting the three features of the complex wavelet result

figure(10), clf

% plot the filtered signal (projection onto real axis)
subplot(311)
plot(timevec,real(convres))
xlabel('Time (ms)'), ylabel('Amplitude (\muV)')
set(gca,'xlim',[-.1 1.4])


% plot power (squared magnitude from origin to dot-product location in complex space)
subplot(312)
plot(timevec,(2*abs(convres)).^2);
xlabel('Time (ms)'), ylabel('Power \muV^2')
set(gca,'xlim',[-.1 1.4])


% plot phase (angle of vector to dot-product, relative to positive real axis)
subplot(313)
plot(timevec,angle(convres))
xlabel('Time (ms)'), ylabel('Phase (rad.)')
set(gca,'xlim',[-.1 1.4])

%% viewing the results as a movie

figure(11), clf

% setup the time course plot
subplot(212)
h = plot(timevec,abs(convres),'k','linew',2);
set(gca,'xlim',timevec([1 end]),'ylim',[min(abs(convres)) max(abs(convres)) ])
xlabel('Time (sec.)'), ylabel('Amplitude (\muV)')


for ti=1:5:length(timevec)
    
    % draw complex values in polar space
    subplot(211)
    polar(0,max(abs(convres))), hold on
    polar(angle(convres(max(1,ti-100):ti)),abs(convres(max(1,ti-100):ti)),'k')
    text(-.75,0,[ num2str(round(1000*timevec(ti))) ' s' ]), hold off
    
    % now show in 'linear' plot
    set(h,'XData',timevec(max(1,ti-100):ti),'YData',abs(convres(max(1,ti-100):ti)))
    
    drawnow
    pause(.00001)
end

%%  VIDEO: Convolution with all trials!
%%% This is the "super-trial" concatenation trick you saw in the slides.
%   Notice the size of the reshaped data and the new convolution parameters.
% extract more data
data = squeeze( csd(6,:,:) );

% reshape the data to be 1D
dataR = reshape(data,1,[]);

% Step 1: N's of convolution
ndata = length(dataR); % note the different variable name!
nkern = length(time);
nConv = ndata + nkern - 1;
halfK = floor(nkern/2);

% Step 2: FFTs
dataX = fft( dataR,nConv );
kernX = fft( cmw,  nConv );

% Step 2.5: normalize the wavelet (try it without this step!)
kernX = kernX ./ max(kernX);

% Step 3: multiply spectra
convresX = dataX .* kernX;

% Step 4: IFFT
convres = ifft(convresX);

% Step 5: cut off "wings"
convres = convres(halfK+1:end-halfK+1);

% New step 6: reshape!
convres2D = reshape(convres,size(data));


%%% now plotting
figure(12), clf
subplot(121)
imagesc(timevec,[],data')
xlabel('Time (s)'), ylabel('Trials')
set(gca,'xlim',[-.1 1.4],'clim',[-1 1]*2000)
title('Broadband signal')

subplot(122)
imagesc(timevec,[],abs(convres2D'))
xlabel('Time (s)'), ylabel('Trials')
set(gca,'xlim',[-.1 1.4],'clim',[-1 1]*500)
title([ 'Power time series at ' num2str(frex) ' Hz' ])


%% VIDEO: A full time-frequency power plot!
%%% Take a deep breath: You're about to make your first time-frequency
%   power plot!

% frequency parameters
min_freq =  5; % in Hz
max_freq = 90; % in HZ
num_freq = 30; % in count

frex = linspace(min_freq,max_freq,num_freq);

% initialize TF matrix
tf = zeros(num_freq,length(timevec));

% IMPORTANT! I'm omitting a few steps of convolution 
%            that are already computed above.

for fi=1:num_freq
    
    % create wavelet
    cmw  = exp(1i*2*pi*frex(fi)*time) .* ...
           exp( -4*log(2)*time.^2 / .3^2 );
    
    cmwX = fft(cmw,nConv);
    cmwX = cmwX./max(cmwX);
    
    % the rest of convolution
    as = ifft( dataX.*cmwX );
    as = as(halfK+1:end-halfK+1);
    as = reshape(as,size(data));
    
    % extract power
    aspow = abs(as).^2;
    
    % average over trials and put in matrix
    tf(fi,:) = mean(aspow,2);
end

%%% and plot!
figure(13), clf
contourf(timevec,frex,tf,40,'linecolor','none')
set(gca,'clim',[0 1]*10000,'xlim',[-.1 1.4])
xlabel('Time (s)'), ylabel('Frequency (Hz)')

%   VIDEO: Inter-trial phase clustering (ITPC/ITC)
% 


%% ITPC with different variances
%%% The goal here is to develop some visual intuition for the
%   correspondence between ITPC and distributions of phase angles.

% specify parameters
circ_prop = .5; % proportion of the circle to fill
N = 100; % number of "trials"

% generate phase angle distribution
simdata = rand(1,N) * (2*pi) * circ_prop;

% compute ITPC and preferred phase angle
itpc      = abs(mean(exp(1i*simdata),2));
prefAngle = angle(mean(exp(1i*simdata)));


% and plot...
figure(14), clf

% as linear histogram
subplot(121)
hist(simdata,20)
xlabel('Phase angle'), ylabel('Count')
set(gca,'xlim',[0 2*pi])
title([ 'Observed ITPC: ' num2str(itpc) ])

% and as polar distribution
subplot(122)
polar([zeros(1,N); simdata],[zeros(1,N); ones(1,N)],'k')
hold on
h = polar([0 prefAngle],[0 itpc],'m');
set(h,'linew',3)
title([ 'Observed ITPC: ' num2str(itpc) ])


%% Compute and plot TF-ITPC for one electrode

load sampleEEGdata.mat

% wavelet parameters
num_frex = 40;
min_freq =  2;
max_freq = 30;

channel2use = 'pz';

% set range for variable number of wavelet cycles
range_cycles = [ 3 10 ];

% parameters (notice using logarithmically spaced frequencies!)
frex  = logspace(log10(min_freq),log10(max_freq),num_frex);
nCycs = logspace(log10(range_cycles(1)),log10(range_cycles(end)),num_frex);
time  = -2:1/EEG.srate:2;
half_wave = (length(time)-1)/2;

% FFT parameters
nWave = length(time);
nData = EEG.pnts*EEG.trials;
nConv = nWave+nData-1;


% FFT of data (doesn't change on frequency iteration)
dataX = fft( reshape(EEG.data(strcmpi(channel2use,{EEG.chanlocs.labels}),:,:),1,nData) ,nConv);

% initialize output time-frequency data
tf = zeros(num_frex,EEG.pnts);

% loop over frequencies
for fi=1:num_frex
    
    % create wavelet and get its FFT
    s = nCycs(fi)/(2*pi*frex(fi));
    wavelet  = exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*s^2));
    waveletX = fft(wavelet,nConv);
    
    % question: is this next line necessary?
    waveletX = waveletX./max(waveletX);
    
    % run convolution
    as = ifft(waveletX.*dataX,nConv);
    as = as(half_wave+1:end-half_wave);
    
    % reshape back to time X trials
    as = squeeze(reshape(as,size(EEG.data(strcmpi(channel2use,{EEG.chanlocs.labels}),:,:))));
    
    % compute ITPC
    tf(fi,:) = abs(mean(exp(1i*angle(as)),2));
end

% plot results
figure(15), clf
contourf(EEG.times,frex,tf,40,'linecolor','none')
set(gca,'clim',[0 .6],'ydir','normal','xlim',[-300 1000])
title('ITPC')
xlabel('Time (ms)'), ylabel('Frequency (Hz)')


%%   VIDEO: Time-frequency trade-off
%%% This code will re-create the plot comparing wavelets with different widths, 
%   in the time and frequency domains.

srate = 512;

% set a few different wavelet widths ("number of cycles" parameter)
num_cycles = [ 2 6 8 15 ];
frex = 8;

time = -2:1/srate:2;
hz = linspace(0,srate/2,floor(length(time)/2)+1);

figure(16), clf
for i=1:4
    
    %%% time domain
    subplot(4,2,i*2-1)
    s = num_cycles(i) / (2*pi* frex); % Gaussian width as number-of-cycles (don't forget to normalize!)
    plot(time, exp( (-time.^2) ./ (2*s^2) ),'k','linew',3)
    title([ 'Gaussian with ' num2str(num_cycles(i)) ' cycles' ])
    xlabel('Time (s)')
    
    
    %%% frequency domain
    subplot(4,2,i*2)
    cmw = exp(1i*2*pi*frex.*time) .* exp( (-time.^2) ./ (2*s^2) );
    
    % take its FFT
    cmwX = fft(cmw);
    cmwX = cmwX./max(cmwX);
    
    % plot it
    plot(hz,abs(cmwX(1:length(hz))),'k','linew',3)
    set(gca,'xlim',[0 20])
    xlabel('Frequency (Hz)')
    title([ 'Power of wavelet with ' num2str(num_cycles(i)) ' cycles' ])
end

%% comparing wavelet convolution with different wavelet settings
clear
load sampleEEGdata.mat

% wavelet parameters
num_frex = 40;
min_freq =  2;
max_freq = 30;

channel2use = 'o1';

% set a few different wavelet widths ("number of cycles" parameter)
num_cycles = [ 2 6 8 15 ];


% time window for baseline normalization.
%  we'll talk about this soon in lecture
baseline_window = [ -500 -200 ];

% other wavelet parameters
frex = linspace(min_freq,max_freq,num_frex);
time = -2:1/EEG.srate:2;
half_wave = (length(time)-1)/2;

% FFT parameters
nKern = length(time);
nData = EEG.pnts*EEG.trials;
nConv = nKern+nData-1;

% initialize output time-frequency data
tf = zeros(length(num_cycles),length(frex),EEG.pnts);

% convert baseline time into indices
baseidx = dsearchn(EEG.times',baseline_window');


% FFT of data (doesn't change on frequency iteration)
dataX = fft( reshape(EEG.data(strcmpi(channel2use,{EEG.chanlocs.labels}),:,:),1,[]) ,nConv);

% loop over cycles
for cyclei=1:length(num_cycles)
    
    for fi=1:length(frex)
        
        % create wavelet and get its FFT
        s = num_cycles(cyclei) / (2*pi*frex(fi));
        
        cmw  = exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*s^2));
        cmwX = fft(cmw,nConv);
        cmwX = cmwX./max(cmwX);
        
        % run convolution, trim edges, and reshape to 2D (time X trials)
        as = ifft(cmwX.*dataX);
        as = as(half_wave+1:end-half_wave);
        as = reshape(as,EEG.pnts,EEG.trials);
        
        % put power data into big matrix
        tf(cyclei,fi,:) = mean(abs(as).^2,2);
    end
    
    % db normalization
    tf(cyclei,:,:) = 10*log10( bsxfun(@rdivide, squeeze(tf(cyclei,:,:)), mean(tf(cyclei,:,baseidx(1):baseidx(2)),3)' ) );
    
end

% plot results
figure(17), clf
for cyclei=1:length(num_cycles)
    subplot(2,2,cyclei)
    
    contourf(EEG.times,frex,squeeze(tf(cyclei,:,:)),40,'linecolor','none')
    set(gca,'clim',[-3 3],'ydir','normal','xlim',[-300 1000])
    title([ 'Wavelet with ' num2str(num_cycles(cyclei)) ' cycles' ])
    xlabel('Time (ms)'), ylabel('Frequency (Hz)')
end

%% variable number of wavelet cycles

% set a few different wavelet widths (number of wavelet cycles)
range_cycles = [ 4 13 ];

% other wavelet parameters
nCycles = logspace(log10(range_cycles(1)),log10(range_cycles(end)),num_frex);

% initialize output time-frequency data
tf = zeros(length(frex),EEG.pnts);

for fi=1:length(frex)
    
    % create wavelet and get its FFT
    s = nCycles(fi)/(2*pi*frex(fi));
    cmw = exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*s^2));
    cmwX = fft(cmw, nConv);
    cmwX = cmwX ./ max(cmwX);
   
    % run convolution
    as = ifft(cmwX .* dataX);
    as = as(half_wave+1:end-half_wave);
    as = reshape(as,EEG.pnts,EEG.trials);
    
    % put power data into big matrix
    tf(fi,:) = mean(abs(as).^2,2);
end

% db normalization (we'll talk about this in the next lecture)
tfDB = 10*log10( bsxfun(@rdivide, tf, mean(tf(:,baseidx(1):baseidx(2)),2)) );

% plot results
figure(18), clf

subplot(2,1,1)
contourf(EEG.times,frex,tf,40,'linecolor','none')
set(gca,'clim',[0 5],'ydir','normal','xlim',[-300 1000])
title('Convolution with a range of cycles')
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
colorbar

subplot(2,1,2)
contourf(EEG.times,frex,tfDB,40,'linecolor','none')
set(gca,'clim',[-3 3],'ydir','normal','xlim',[-300 1000])
title('Same data but dB normalized!')
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
colorbar
%%   VIDEO: Baseline normalization of TF plots
% specify baseline periods for dB-normalization
baseline_windows = [ -500 -200;
                     -100    0;
                        0  300;
                     -800    0;
                   ];

               
% convert baseline time into indices
baseidx = reshape( dsearchn(EEG.times',baseline_windows(:)), [],2);

%% setup wavelet parameters

% frequency parameters
min_freq =  2;
max_freq = 30;
num_frex = 40;
frex = linspace(min_freq,max_freq,num_frex);

% which channel to plot
channel2use = 'o1';

% other wavelet parameters
range_cycles = [ 4 10 ];

% notice: defining cycles as a vector for all frequencies
s = logspace( log10(range_cycles(1)) , log10(range_cycles(end)), num_frex) ./ (2*pi*frex);
wavtime = -2:1/EEG.srate:2;
half_wave = (length(wavtime)-1)/2;


% FFT parameters
nWave = length(wavtime);
nData = EEG.pnts * EEG.trials;
nConv = nWave + nData - 1;


% now compute the FFT of all trials concatenated
alldata = reshape( EEG.data(strcmpi(channel2use,{EEG.chanlocs.labels}),:,:) ,1,[]);
dataX   = fft( alldata ,nConv );


% initialize output time-frequency data
tf = zeros(size(baseidx,1),length(frex),EEG.pnts);

%% now perform convolution

% loop over frequencies
for fi=1:length(frex)
    
    % create wavelet and get its FFT
    wavelet  = exp(2*1i*pi*frex(fi).*wavtime) .* exp(-wavtime.^2./(2*s(fi)^2));
    waveletX = fft(wavelet,nConv);
    waveletX = waveletX ./ max(waveletX);
    
    % now run convolution in one step
    as = ifft(waveletX .* dataX);
    as = as(half_wave+1:end-half_wave);
    
    % and reshape back to time X trials
    as = reshape( as, EEG.pnts, EEG.trials );
    
    % compute power and average over trials
    tf(4,fi,:) = mean( abs(as).^2 ,2);
end

%% db normalization and plot results

% define color limits
clim = [-3 3];

% create new matrix for percent change
tfpct = zeros(size(tf));

for basei=1:size(tf,1)
    
    activity = tf(4,:,:);
    baseline = mean( tf(4,:,baseidx(basei,1):baseidx(basei,2)) ,3);
    
    % decibel
    tf(basei,:,:) = 10*log10( activity ./ baseline);
end


% plot results
figure(19), clf
for basei=1:size(baseline_windows,1)
    
    subplot(2,2,basei)
    
    contourf(EEG.times,frex,squeeze(tf(basei,:,:)),40,'linecolor','none')
    set(gca,'clim',clim,'ydir','normal','xlim',[-300 1000])
    title([ 'DB baseline of ' num2str(baseline_windows(basei,1)) ' to ' num2str(baseline_windows(basei,2)) ' ms' ])
end

xlabel('Time (ms)'), ylabel('Frequency (Hz)')


