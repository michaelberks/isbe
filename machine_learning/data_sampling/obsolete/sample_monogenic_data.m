function [mono_samples] = sample_monogenic_data(local_amp, local_phase, local_ori, rows, cols, varargin)
%SAMPLE_MONOGENIC_DATA *Insert a one line summary here*
%   [dt_samples] = sample_dt_polar13(dt, rows, cols, rotate)
%
% Inputs:
%      image_in - image to decompose and sample
%
%      num_levels - num_levels to include
%
%      rows - row subscripts of sample locations
%
%      cols - column subscripts of sample locations
%
%      win_size - dimensions of window to sample
%
%      rotate - if true (default false) circular shift the columns of the
%      final representation so maximal orientation is in first row, thus
%      achieving approximate rotation invariance (+/- 15 degrees)
%
%
% Outputs:
%      linop_samples - *Insert description of input variable here*
%
%
% Example:
%
% Notes: Representation is based on that presented by Kingsbury in
% "Rotation-Invariant Local Feature Matching with Complex Wavelets"
%
% See also:
%
% Created: 23-Aug-2010
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester

% Unpack the arguments:
args = u_packargs(varargin, '0',...
    'win_size', 3,...
    'pca', []);

win_size = args.win_size; 
num_levels = size(local_amp,3)-1;

%Make copies of sample rows and cols at positions of local window patch
pad_w = floor(win_size/2);
win_idx = -pad_w:pad_w;
num_samples = length(rows);
rr = repmat(rows*ones(1,win_size) + ones(num_samples,1)*win_idx, 1, win_size);
cc = kron(cols*ones(1,win_size) + ones(num_samples,1)*win_idx, ones(1,win_size));

patch_idx = sub2ind(size(local_amp(:,:,1)), rr, cc);

%save sample into output container
local_amp = reshape(local_amp(:,:,2:end), [], num_levels);
local_phase = reshape(local_phase(:,:,2:end), [], num_levels);
local_ori = reshape(local_ori(:,:,2:end), [], num_levels);

mono_samples = zeros(num_samples, 3*num_levels*win_size^2);

for level = 1:num_levels
    amp_cols = (1:win_size^2) + (level-1)*win_size^2;
    mono_samples(:, amp_cols) = ...
        reshape(local_amp(patch_idx,level), num_samples, win_size^2);

    phase_cols = (1:win_size^2) + (num_levels + level-1)*win_size^2;
    mono_samples(:, phase_cols) = ...
        reshape(local_phase(patch_idx,level),num_samples, win_size^2);
    
    ori_cols = (1:win_size^2) + (2*num_levels + level-1)*win_size^2;
    mono_samples(:, ori_cols) = ...
        reshape(local_ori(patch_idx,level),num_samples, win_size^2);
end

if ~isempty(args.pca)
    %Transform sample using PCA modes
    mono_samples = bsxfun(@minus, mono_samples, args.pca.mean)*args.pca.modes;
end