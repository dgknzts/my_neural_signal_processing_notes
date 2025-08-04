%% CSD Analysis Tutorial: Using CSD Toolbox
% Demonstrates spatial sharpening with professional CSD implementation
%
% Uses CSD Toolbox for MATLAB:
% Kayser, J., & Tenke, C. E. (2006). Principal components analysis of
% Laplacian waveforms as a generic method for identifying ERP generator
% patterns: I. Evaluation with auditory oddball tasks. Clinical
% Neurophysiology, 117(2), 348-368.
%
% Toolbox available at: http://psychophysiology.cpmc.columbia.edu/software/csdtoolbox/
clear; close all;
addpath(genpath('your path to the CSD toolbox')); % CSD toolbox path

%% 1. CRITICAL DATA IMPORT CONSIDERATIONS
% Before applying CSD, ensure:
% - Accurate 3D electrode coordinates matching CSD toolbox files
% - Remove non-scalp channels (EOG, EMG, reference electrodes)
% - Artifact-free data (CSD amplifies noise via spatial differentiation)
% - Minimum 31 electrodes (64+ preferred for accurate spatial sampling)
% - Proper filtering and baseline correction applied
% - Screen for electrolyte bridges between adjacent electrodes
% - Any reference scheme acceptable (CSD produces reference-free output)

%% 2. Import your data with chan x time x epochs

%% 3. APPLY CSD TRANSFORMATION USING TOOLBOX
% Extract electrode labels for CSD toolbox
electrode_labels = {EEG.chanlocs.labels}';

% Remove non-scalp channels if any
scalp_idx = ~strcmpi(electrode_labels, 'VEOG') & ...
 ~strcmpi(electrode_labels, 'HEOG') & ...
 ~contains(electrode_labels, 'EMG');
if sum(~scalp_idx) > 0
 EEG = pop_select(EEG, 'channel', find(scalp_idx));
 electrode_labels = {EEG.chanlocs.labels}';
 n_channels = EEG.nbchan;
end

% Extract montage using CSD toolbox
% Matches electrode labels to predefined 3D coordinates
[M] = ExtractMontage('10-5-System_Mastoids_EGI129.csd', electrode_labels);

% CSD parameters
m_constant = 4; % Spline flexibility (4 = recommended balance)
lambda = 1e-5; % Smoothing regularization parameter
head_radius = 10; % Head radius for scaling (cm)

% Generate CSD transformation matrices
% G: Cosine matrix of inter-electrode angular distances on sphere
% H: Spherical spline coefficients from Legendre polynomials
[G, H] = GetGH(M, m_constant);

% Apply CSD transformation to each epoch
EEG_csd = EEG; % Copy structure
csd_data = zeros(size(EEG.data));

for epoch = 1:n_epochs
% CSD function computes negative second spatial derivative
% Estimates radial current flow: sources (+) toward scalp, sinks (-) away
% Output: reference-free current density (μV/cm²)
 epoch_data = squeeze(EEG.data(:, :, epoch));
 csd_data(:, :, epoch) = CSD(epoch_data, G, H, lambda, head_radius);
end

EEG_csd.data = csd_data;