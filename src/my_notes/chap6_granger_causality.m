%%  VIDEO: Granger causality (or prediction, whatever)
% univariate autoregressive model parameters
a1 = 1.2;
a2 = -.9;
N = 30;

% start with two random points
x = randn(2,1);

% add new points according to previous weighted values
for t=2:N-1
    x(t+1) = a1*x(t) + a2*x(t-1);
end


% and plot
figure(11), clf
plot(x,'s-','markerfacecolor','k','markersize',10,'linew',2)


%% bivariate autoregression

% regression parameters
a1 = 1.2;
a2 = -.9;
N = 30;

% initialize data vectors
x = randn(N,1);
y = randn(2,1);

% new y
for t=2:N-1
    y(t+1) = a1*x(t) + a2*x(t-1);
end

figure(12), clf, hold on
plot(1:N,x  ,'s-','markerfacecolor','k','markersize',10,'linew',2)
plot(1:N,y-3,'s-','markerfacecolor','k','markersize',10,'linew',2)
legend({'x';'y'},'box','off')
set(gca,'ytick',[])

%% a longer example of bivariate autoregression

% frequencies for the sine waves
freq1 = 10;
freq2 = 2;

% simulation parameters
srate = 1000;
time  = 0:1/srate:1;

% define x
x = .39*sin(2*pi*freq1.*time) + .7*sin(2*pi*freq2.*time) + randn(size(time))/5;

% define y by x
y = [0 .2];
for i=3:length(x)
    y(i) = .8*x(i-2) + 1.2*x(i-1) + randn/5;
end

figure(13), clf
plot(time,x,'m',time,y,'k')

xlabel('Time (ms)')
ylabel('Amplitude (arb. units)')
title('Bivariate autoregression')

legend({[ 'x = .39*sin(2\pi' num2str(freq1) 't) + .7*sin(2\pi' num2str(freq2) 't) + \sigma' ],'y = .8x_{t-2} + 1.2x_{t-1} + \sigma'})


%% Autoregression model estimation
% construct autoregression model
x=1;
for i=1:30
    x(i+1) = 1.1*x(i);
end

% recover parameters using armorf function (note: this function is taken
% from the BSMART toolbox (available online).
[Ax,Ex] = armorf(x,1,length(x),1);
fprintf('\n\n')
disp('Univariate AR, timeseries 1:')
disp([ '  Autoregression coefficient (order=1): ' num2str(Ax) ', Error: ' num2str(Ex) ]);

% try using order=2;
[Ax,Ex] = armorf(x,1,length(x), 2);
disp([ '  Autoregression coefficients (order=2): ' num2str(Ax) ', Error: ' num2str(Ex) ]);
disp(' ')





%%% another example with true order=2
x=[1 1.5];
for i=1:30
    x(i+2) = 1.1*x(i) + -.3*x(i+1);
end

% recover parameters using armorf function
[Ax,Ex] = armorf(x,1,length(x),1);
disp('Univariate AR, timeseries 2:')
disp([ '  Autoregression coefficient (order=1): ' num2str(Ax) ', Error: ' num2str(Ex) ]);

% try using order=2;
[Ax,Ex] = armorf(x,1,length(x),2);
disp([ '  Autoregression coefficients (order=2): ' num2str(Ax) ', Error: ' num2str(Ex) ]);
fprintf('\n')






%%% another example with bivariate
x = .39*sin(2*pi*freq1.*time) + .7*sin(2*pi*freq2.*time) + randn(size(time))/10;
y = [0 0];
for i=1:length(x)-2
    y(i+2) = -.8*x(i) + 1.2*x(i+1);
end

% 
[Axy,Exy] = armorf([x; y],1,length(x),2);
disp('Bivariate AR:')
Axy
Exy

%%%%% NOTE: In the video I explained the organization of Axy incorrectly.
%%%%% The correct organization is as follows:
% Note: variable Axy is organized as:
%    XX_1 XX_1  YX_2 YX_2
%    XY_1 XY_1  YY_2 YY_2


%% Directed synchrony through Granger prediction


% load our favorite EEG data
load sampleEEGdata.mat

% define channels to compute granger synchrony between
chan1name = 'fcz';
chan2name = 'o1';

% find the index of those channels
chan1 = find( strcmpi(chan1name,{EEG.chanlocs.labels}) );
chan2 = find( strcmpi(chan2name,{EEG.chanlocs.labels}) );

% define autoregression parameters (can leave as default for now)
order = 14;


% get AR coefficients and error from each signal
[Ax,Ex] = armorf(EEG.data(chan1,:,1),1,EEG.pnts,order);
[Ay,Ey] = armorf(EEG.data(chan2,:,1),1,EEG.pnts,order);


%%% we are going to reconstruct the data using the autoregressive coefficients
x = zeros(1,EEG.pnts);
y = zeros(1,EEG.pnts);

x(1:order) = EEG.data(chan1,1:order,1);
y(1:order) = EEG.data(chan2,1:order,1);

for i=order+1:EEG.pnts
    
    % initialize
    thispointX = 0;
    thispointY = 0;
    
    for ai=1:order
        thispointX = thispointX + EEG.data(chan1,i-ai,1)*Ax(ai);
        thispointY = thispointY + EEG.data(chan2,i-ai,1)*Ay(ai);
    end
    x(i-1) = thispointX;
    y(i-1) = thispointY;
end

figure(14), clf
subplot(211)
plot(EEG.times,EEG.data(chan1,:,1),'b', EEG.times,x,'r')
legend({'Real data';'Reconstructed from ARmodel'})

subplot(212)
plot(EEG.times,EEG.data(chan2,:,1),'b', EEG.times,y,'r')
legend({'Real data';'Reconstructed from ARmodel'})

%% Now for Granger prediction


% Bivariate autoregression and associated error term
[Axy,E] = armorf(EEG.data([chan1 chan2],:,1),1,EEG.pnts,order);


% time-domain causal estimate
granger_chan2_to_chan1 = log(Ex/E(1,1));
granger_chan1_to_chan2 = log(Ey/E(2,2));

disp([ 'Granger prediction from ' chan1name ' to ' chan2name ' is ' num2str(granger_chan1_to_chan2) ]);
disp([ 'Granger prediction from ' chan2name ' to ' chan1name ' is ' num2str(granger_chan2_to_chan1) ]);

%% Now we compute granger prediction over time

% initialize
x2yT = zeros(1,EEG.pnts);
y2xT = zeros(1,EEG.pnts);

% GC parameters
iwin   = 300; % in ms
iorder = 15;  % in ms


% convert window/order to points
win   = round(iwin/(1000/EEG.srate));
order = round(iorder/(1000/EEG.srate));

for timei=1:EEG.pnts-win
    
    % data from all trials in this time window
    % Data should be normalized before computing Granger estimates
    tempdata = zscore(reshape(EEG.data([chan1 chan2],timei:timei+win-1,1),2,win),0,2);
    
    %% fit AR models (model estimation from bsmart toolbox)
    [Ax,Ex] = armorf(tempdata(1,:),1,win,order);
    [Ay,Ey] = armorf(tempdata(2,:),1,win,order);
    [Axy,E] = armorf(tempdata     ,1,win,order);
    
    % time-domain causal estimate
    y2xT(timei) = log(Ex/E(1,1));
    x2yT(timei) = log(Ey/E(2,2));
    
end

% draw lines
figure(15), clf, hold on

plot(EEG.times,x2yT)
plot(EEG.times,y2xT,'r')
legend({[ 'GC: ' chan1name ' -> ' chan2name ];[ 'GC: ' chan2name ' -> ' chan1name ]})

title([ 'Window length: ' num2str(iwin) ' ms, order: ' num2str(iorder) ' ms' ])
xlabel('Time (ms)')
ylabel('Granger prediction estimate')
set(gca,'xlim',[-200 1000])
