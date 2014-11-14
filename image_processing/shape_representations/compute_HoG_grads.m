function [g_mag g_theta_idx] = compute_HoG_grads(image_in, varargin)
%COMPUTE_HOG *Insert a one line summary here*
%   [] = compute_HoG(varargin)
%
% COMPUTE_HOG uses the U_PACKARGS interface function
% and so arguments to the function may be passed as name-value pairs
% or as a struct with fields with corresponding names. The field names
% are defined as:
%
% Mandatory Arguments:
%
% Optional Arguments:
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 30-Sep-2013
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

% Unpack the arguments:
args = u_packargs(varargin, '0', ...
    'cell_sz', [8 8],... %Size of cells in blocks
    'block_sz', [2 2],...%Size of blocks in cells
    'num_ori_bins', 9,... %Number of bins in orientation histograms
    'norm_method', 'l1-sqrt',... %Method for local normalisation
    'block_spacing', 8,... %Separation between blocks in pixels - controls level of overlap
    'gradient_operator', [-1 0 1],...
    'spatial_sigma', 0, ...
    'window_sz', [64 64],...
    'angle_wrap', 0,...
    'debug', 0);
clear varargin;

%Apply gradients to input image
gy = conv2(image_in, args.gradient_operator', 'same');
gx = conv2(image_in, args.gradient_operator, 'same');

g_mag = sqrt(gx.^2 + gy.^2);

if args. angle_wrap
    %Compute angle on half-circle (will be [-pi/2 pi/2]) then scale to [0
    %1]
    g_theta_idx = atan(g_y ./ g_x)/pi + 0.5;
else
    %Compute angle on full circle (will be [-pi pi]) then scale to [0 1]
    g_theta_idx = atan2(g_y, g_x) / (2*pi) + 0.5;
end

%Convert theta vals in into integer indices - there may still be some 0s,
%which should be moved to the final bin
g_theta_idx = ceil(g_theta_idx * args.num_ori_bins);
g_theta_idx(~g_theta_idx) = args.num_ori_bins;


