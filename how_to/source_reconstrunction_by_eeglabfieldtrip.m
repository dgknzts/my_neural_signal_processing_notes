%% EEG Source Localization Tutorial using FieldTrip
% This tutorial demonstrates how to perform source localization on EEG data
% using the DICS (Dynamic Imaging of Coherent Sources) beamformer method
%
% Main steps:
% 1. Load and prepare EEG data
% 2. Preprocess data (rereferencing)
% 3. Define pre/post stimulus periods
% 4. Frequency analysis using wavelets
% 5. Forward modeling (leadfield calculation)
% 6. Source analysis using DICS beamformer
% 7. Visualization of results

%% 1. SETUP AND DATA LOADING
% Initialize toolboxes
run("G:\My Drive\Projects\toolboxes\eeglab2025.0.0\eeglab.m")
addpath("G:\My Drive\Projects\toolboxes\fieldtrip-20250106\fieldtrip-20250106")
ft_defaults

%% Load EEG dataset with simulated dipole activity at 15Hz
EEG = pop_loadset('filename','ICAsimulated_EEG_freq15Hz_phase0_dip56_32.set',...
                  'filepath','G:\\My Drive\\Projects\\signal_processing_mike_cohen\\');

% Set electrode positions using standard 10-05 system
EEG = pop_chanedit(EEG, 'lookup',...
    'G:\\My Drive\\Projects\\toolboxes\\eeglab2025.0.0\\plugins\\dipfit\\standard_BEM\\elec\\standard_1005.elc');

% Configure dipole fitting parameters for source localization
EEG = pop_dipfit_settings(EEG, 'hdmfile','standard_vol.mat',...
                         'mrifile','standard_mri.mat',...
                         'chanfile','standard_1005.elc',...
                         'coordformat','MNI',...
                         'coord_transform',[-1.8347e-06 -9.9255e-06 1.4786e-05 1.099e-07 -1.2661e-07 -1.5708 1 1 1]);

% Convert EEGLAB structure to FieldTrip format
data_all = eeglab2fieldtrip(EEG, 'preprocessing', 'dipfit');


%% 2. PREPROCESSING
% Apply common average reference to reduce volume conduction effects
cfg = [];
cfg.channel = {'all'};
cfg.reref = 'yes';           % Enable rereferencing
cfg.refchannel = {'all'};    % Use all channels for common average reference
data_all = ft_preprocessing(cfg, data_all);

fprintf('Applied common average reference to reduce volume conduction\n');

%% 3. DEFINE PRE/POST STIMULUS PERIODS
% Split data into baseline (pre-stimulus) and activation (post-stimulus) periods
% This allows comparison of brain activity before and after stimulation

% Pre-stimulus period: 0.5 to 0.8 seconds (baseline)
cfg = [];
cfg.latency = [0.5 0.8-1./data_all.fsample];  % Subtract one sample to avoid overlap
dataPre = ft_selectdata(cfg, data_all);

% Post-stimulus period: 1.0 to 3.0 seconds (activation)
cfg.latency = [1 3-1./data_all.fsample];
dataPost = ft_selectdata(cfg, data_all);

fprintf('Data segmented: Pre-stim (%.1f-%.1fs), Post-stim (%.1f-%.1fs)\n', ...
    dataPre.time{1}(1), dataPre.time{1}(end), ...
    dataPost.time{1}(1), dataPost.time{1}(end));

%% 4. FREQUENCY ANALYSIS
% Compute power spectral density and cross-spectral density for beamforming

% Post-stimulus: Use wavelet analysis for time-frequency resolution
cfg = [];
cfg.method = 'wavelet';      % Wavelet analysis
cfg.width = 7;              % Wavelet width (cycles)
cfg.output = 'powandcsd';    % Output both power and cross-spectral density
cfg.foi = 15;               % Frequency of interest (15 Hz)
cfg.toi = 1.5:0.5:2.5;     % Time points of interest (every 50ms)
freqPost = ft_freqanalysis(cfg, dataPost);

fprintf('Frequency analysis completed for 15 Hz signal\n');

%% 5. FORWARD MODEL SETUP
% Load standard head model for source localization

% Load boundary element method (BEM) head model
headmodel = ft_read_headmodel('standard_bem.mat');

% Visualize head model (brain, skull, scalp boundaries)
figure('Name', 'Head Model');
ft_plot_headmodel(headmodel, 'facealpha', 0.6);
title('Standard BEM Head Model');
camlight; lighting gouraud;

% Load standard MRI template
mri = ft_read_mri('standard_mri.mat');

% Visualize MRI with head model intersection
figure('Name', 'MRI with Head Model');
cfg = [];
cfg.intersectmesh = {headmodel.bnd(1), headmodel.bnd(2), headmodel.bnd(3)};
ft_sourceplot(cfg, mri);
title('MRI Template with Head Model Boundaries');

fprintf('Head model and MRI template loaded\n');

%% 6. SOURCE MODEL PREPARATION
% Create 3D grid for source reconstruction

cfg = [];
cfg.elec = freqPost.elec;    % Use electrode positions from data
cfg.headmodel = headmodel;   % Use BEM head model
cfg.reducerank = 2;          % Reduce rank for EEG (no radial component)
cfg.resolution = 1;          % 1 cm resolution grid
cfg.sourcemodel.unit = 'cm'; % Units in centimeters
cfg.normalize = 'yes';       % Normalize leadfield
sourcemodel = ft_prepare_leadfield(cfg);

fprintf('Leadfield computed: %d sources in brain volume\n', ...
    sum(sourcemodel.inside));

%% 7. SOURCE ANALYSIS USING DICS BEAMFORMER
% Apply Dynamic Imaging of Coherent Sources (DICS) beamformer

cfg = [];
cfg.method = 'dics';         % DICS beamformer
cfg.frequency = 15;          % Frequency of interest
cfg.sourcemodel = sourcemodel;
cfg.headmodel = headmodel;
cfg.dics.projectnoise = 'yes';  % Project noise for normalization
cfg.dics.lambda = 0;         % No regularization (can be increased if unstable)
sourcePost_nocon = ft_sourceanalysis(cfg, freqPost);

fprintf('DICS beamformer analysis completed\n');

%% 8. SOURCE INTERPOLATION
% Interpolate source results onto MRI for visualization

% Reslice MRI to match coordinate system
cfg = [];
mri = ft_volumereslice(cfg, mri);

% Interpolate source power onto MRI
cfg = [];
cfg.downsample = 2;          % Downsample for faster visualization
cfg.parameter = 'pow';       % Interpolate power values
sourcePostInt_nocon = ft_sourceinterpolate(cfg, sourcePost_nocon, mri);

%% 9. VISUALIZATION OF RAW POWER
% Display source power without normalization

figure('Name', 'Raw Source Power');
cfg = [];
cfg.method = 'slice';        % Slice-based visualization
cfg.funparameter = 'pow';    % Display power
ft_sourceplot(cfg, sourcePostInt_nocon);
title('Raw Source Power at 15 Hz');

%% 10. NEURAL ACTIVITY INDEX (NAI) CALCULATION
% Normalize power by noise to enhance signal-to-noise ratio

sourceNAI = sourcePost_nocon;
sourceNAI.avg.pow = sourcePost_nocon.avg.pow ./ sourcePost_nocon.avg.noise;

% Interpolate NAI onto MRI
cfg = [];
cfg.downsample = 2;
cfg.parameter = 'pow';       % Now contains NAI values
sourceNAIInt = ft_sourceinterpolate(cfg, sourceNAI, mri);

fprintf('Neural Activity Index (NAI) calculated\n');

%% 11. FINAL VISUALIZATION WITH THRESHOLDING
% Display NAI with appropriate thresholds and transparency

maxval = max(sourceNAIInt.pow, [], 'all');

figure('Name', 'Neural Activity Index (NAI)');
cfg = [];
cfg.method = 'slice';
cfg.funparameter = 'pow';         % Display NAI
cfg.maskparameter = cfg.funparameter;  % Use same parameter for masking
cfg.funcolorlim = [2.6 maxval];   % Color limits (threshold at NAI=3)
cfg.opacitylim = [2.6 maxval];    % Opacity limits (transparent below NAI=3)
cfg.opacitymap = 'rampup';        % Gradual opacity increase
ft_sourceplot(cfg, sourceNAIInt);
title(sprintf('Neural Activity Index (Threshold: %.1f, Max: %.1f)', 2, maxval));

%% ADDITIONAL ANALYSIS SUGGESTIONS
% 
% To extend this analysis, consider:
% 1. Compare multiple frequency bands
% 2. Perform statistical testing across trials
% 3. Apply different beamformer variants (LCMV, SAM)
% 4. Investigate connectivity between sources
% 5. Validate results with known simulated dipole locations