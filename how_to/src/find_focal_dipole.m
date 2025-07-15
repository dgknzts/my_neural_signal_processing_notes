function [dipole_idx, focality_score, sorted_dipoles] = find_focal_dipole(lf, chanlocs, target_channels, varargin)
% Find dipoles with focal projection to specific channels
%
% Inputs:
%   lf              - leadfield structure with .Gain field
%   chanlocs        - EEG channel locations structure
%   target_channels - string or cell array of target channel names
%   'method'        - 'ratio' (default) or 'zscore'
%   'n_neighbors'   - number of nearby channels to include as targets (default: 0)
%
% Outputs:
%   dipole_idx      - index of most focal dipole
%   focality_score  - focality score
%   sorted_dipoles  - all dipoles sorted by focality

% Parse inputs
p = inputParser;
addParameter(p, 'method', 'ratio', @ischar);
addParameter(p, 'n_neighbors', 0, @isnumeric);
parse(p, varargin{:});
method = p.Results.method;
n_neighbors = p.Results.n_neighbors;

% Handle single channel input
if ischar(target_channels)
    target_channels = {target_channels};
end

% Find target channel indices
n_targets = length(target_channels);
target_idx = zeros(n_targets, 1);
for i = 1:n_targets
    idx = find(strcmpi({chanlocs.labels}, target_channels{i}));
    if isempty(idx)
        error('Channel %s not found', target_channels{i});
    end
    target_idx(i) = idx;
end

% Add neighboring channels if requested
if n_neighbors > 0
    target_idx = add_neighbors(target_idx, chanlocs, n_neighbors);
end

% Non-target indices
all_idx = 1:length(chanlocs);
nontarget_idx = setdiff(all_idx, target_idx);

% Calculate focality for each dipole
n_dipoles = size(lf.Gain, 3);
focality = zeros(n_dipoles, 1);

for d = 1:n_dipoles
    % Projection strength to target channels
    target_gain = squeeze(lf.Gain(target_idx, :, d));
    if numel(target_idx) == 1
        target_gain = target_gain(:)';
    end
    target_strength = sqrt(mean(target_gain(:).^2));
    
    % Projection strength to non-target channels
    nontarget_gain = squeeze(lf.Gain(nontarget_idx, :, d));
    nontarget_strength = sqrt(mean(nontarget_gain(:).^2));
    
    % Calculate focality score
    switch method
        case 'ratio'
            % Higher ratio = more focal
            focality(d) = target_strength / (nontarget_strength + eps);
        case 'zscore'
            % Z-score of target relative to all channels
            all_strengths = sqrt(squeeze(mean(lf.Gain(:,:,d).^2, 2)));
            focality(d) = (target_strength - mean(all_strengths)) / std(all_strengths);
    end
end

% Find most focal dipole
[focality_score, dipole_idx] = max(focality);

% Sort all dipoles by focality
[sorted_scores, sorted_idx] = sort(focality, 'descend');
sorted_dipoles.indices = sorted_idx;
sorted_dipoles.scores = sorted_scores;
sorted_dipoles.target_idx = target_idx;
sorted_dipoles.method = method;

end

function expanded_idx = add_neighbors(target_idx, chanlocs, n_neighbors)
% Add spatially neighboring channels to target indices

% Get 3D positions
pos = [[chanlocs.X]', [chanlocs.Y]', [chanlocs.Z]'];

expanded_idx = target_idx;
for t = 1:length(target_idx)
    % Calculate distances from target channel
    target_pos = pos(target_idx(t), :);
    distances = sqrt(sum((pos - target_pos).^2, 2));
    
    % Find nearest neighbors
    [~, sorted_idx] = sort(distances);
    neighbors = sorted_idx(2:n_neighbors+1); % exclude self
    
    expanded_idx = unique([expanded_idx; neighbors]);
end

end