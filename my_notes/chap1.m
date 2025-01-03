%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Signal Processing Notes of Mike Cohen %%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Chapter 1: The Basics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Drawing some plots
% importing the example data
load("data/chap1/sampleEEGdata.mat");

% to reach the functions
addpath(genpath("data"));
% Drawing ERP
%% plot ERPs and topographical maps

% compute the ERP of each channel 
erp = mean(EEG.data, 3);


% pick a channel and plot ERP
chan2plot = 'fcz';

figure(1), clf
plot(EEG.times,erp( strcmpi({EEG.chanlocs.labels},chan2plot) ,:),'linew',2)
xlabel('Time (ms)'), ylabel('Activity (\muV)')
set(gca,'xlim',[-400 1200])
%Drawing topos
time2plot = 500; % in ms

% convert time in ms to time in indices
[~,tidx] = min(abs(EEG.times-time2plot));
%min function returns a minimum value and the index of that value. we only
%saving index here !!

% make a topographical map of the ERP at only this time point.
figure(2), clf
topoplotIndie(erp(:,tidx),EEG.chanlocs);
title([ 'ERP from ' num2str(time2plot) ' ms' ])
colorbar
set(gca, 'clim', [-8 8])

%% Other dataset
%Importing other example dataset
load('data/chap1/v1_laminar.mat')
%Multiplechannels heatmap drawing. first 16 channels
% plot depth-by-time image of ERP by averaging over trials
figure(3), clf
contourf(timevec,1:16,squeeze(mean(csd,3)),40,'linecolor','none')
set(gca,'xlim',[0 1.3])
xlabel('Time (sec.)'), ylabel('Cortical depth')

%% Simulating Data
% Check the code for example simulations 
% data/chap1/uANTS_simulate_problemset_SOL.m