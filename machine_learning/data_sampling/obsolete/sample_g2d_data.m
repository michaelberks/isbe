function [g2d_samples] = sample_g2d_data(gf1, gf2, gf3, rows, cols, varargin)
%SAMPLE_G2D_DATA *Insert a one line summary here*
%   [dt_samples] = sample_g2d_data(dt, rows, cols, rotate)
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
%
%
% Example:
%
%
% See also:
%
% Created: 23-Aug-2010
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester

f_debug = (nargin == 0 && nargout == 0);
if f_debug, test_script(); return; end

% Call new function
g2d_samples = sample_filter_responses(cat(4, gf1, gf2, gf3), ...
                             rows, cols, varargin{:});
                         

%% Old code (for posterity)
function g2d_samples = old_code(gf1, gf2, gf3, rows, cols, varargin)

% Unpack the arguments:
args = u_packargs(varargin, '0',...
    'win_size', 3,...
    'pca', [], ...
	'sigma_range', []); % not used but prevents warnings

win_size = args.win_size; 
num_levels = size(gf1,3);

%Make copies of sample rows and cols at positions of local window patch
pad_w = floor(win_size/2);
win_idx = -pad_w:pad_w;
num_samples = length(rows);
rr = repmat(rows*ones(1,win_size) + ones(num_samples,1)*win_idx, 1, win_size);
cc = kron(cols*ones(1,win_size) + ones(num_samples,1)*win_idx, ones(1,win_size));

patch_idx = sub2ind(size(gf1(:,:,1)), rr, cc);

%save sample into output container
gf1 = reshape(gf1, [], num_levels);
gf2 = reshape(gf2, [], num_levels);
gf3 = reshape(gf3, [], num_levels);

g2d_samples = zeros(num_samples, 3*num_levels*win_size^2);

for level = 1:num_levels
    cols1 = (1:win_size^2) + (level-1)*win_size^2;
    g2d_samples(:, cols1) = ...
        reshape(gf1(patch_idx,level), num_samples, win_size^2);

    cols2 = (1:win_size^2) + (num_levels + level-1)*win_size^2;
    g2d_samples(:, cols2) = ...
        reshape(gf2(patch_idx,level),num_samples, win_size^2);
    
    cols3 = (1:win_size^2) + (2*num_levels + level-1)*win_size^2;
    g2d_samples(:, cols3) = ...
        reshape(gf3(patch_idx,level),num_samples, win_size^2);
end

if ~isempty(args.pca)
    %Transform sample using PCA modes
    g2d_samples = bsxfun(@minus, g2d_samples, args.pca.mean)*args.pca.modes;
end


%% Script to test that output of new code matches that of the old code
function test_script()
clc;

% Generate dummy data and compare outputs with sample_filter_responses
n = 5;
rows = 2:n-1; rows = rows(:);
cols = 2:n-1; cols = cols(:);

gf1 = randn(n,n,n);
gf2 = randn(n,n,n);
gf3 = randn(n,n,n);

old_samples = old_code(gf1, gf2, gf3, rows, cols);
new_samples = sample_filter_responses(cat(4,gf1,gf2,gf3), rows, cols);

if any(old_samples ~= new_samples)
    error;
else
    display(['Test successful']);
end
