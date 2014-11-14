function [hog g_mag g_theta_idx] = compute_HoG(image_in, varargin)
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
    'angle_wrap', 1,...
    'debug', 0);
clear varargin;

%Apply gradients to input image
if args.spatial_sigma
    g = gaussian_filters_1d(args.spatial_sigma);
    g = g / sum(g);
    image_in = conv2(g', g, image_in, 'same');
end
g_y = conv2(image_in, args.gradient_operator', 'same');
g_x = conv2(image_in, args.gradient_operator, 'same');

pad_sz = floor(length(args.gradient_operator)/2);
g_x([1:pad_sz end:end-pad_sz+1], :) = [];
g_y([1:pad_sz end:end-pad_sz+1], :) = [];
g_x(:,[1:pad_sz end:end-pad_sz+1]) = [];
g_y(:,[1:pad_sz end:end-pad_sz+1]) = [];

g_mag = sqrt(g_x.^2 + g_y.^2);

if args. angle_wrap
    %Compute angle on half-circle (will be [-pi/2 pi/2]) then scale to [0
    %1]
    g_theta_idx = atan(g_y ./ g_x)/pi + 0.5;
else
    %Compute angle on full circle (will be [-pi pi]) then scale to [0 1]
    g_theta_idx = atan2(g_y, g_x) / (2*pi) + 0.5;
end
g_theta_idx(~g_mag) = 1;

%Convert theta vals in into integer indices - there may still be some 0s,
%which should be moved to the final bin
g_theta_idx = ceil(g_theta_idx * args.num_ori_bins);
g_theta_idx(~g_theta_idx) = args.num_ori_bins;

%Make cell index array
[n_im_rows n_im_cols] = size(g_x);
im_rows = repmat((1:n_im_rows)', 1, n_im_cols);
im_cols = repmat(1:n_im_cols, n_im_rows, 1);

cl_rows = ceil(im_rows / args.cell_sz(1));
cl_cols = ceil(im_cols / args.cell_sz(2));
n_cl_rows = cl_rows(end);
n_cl_cols = cl_cols(end);

cl_idx = sub2ind([n_cl_rows n_cl_cols], cl_rows, cl_cols);

%Now loop through each cell and return the magnitude wieghted orientation
%histogram
num_cells = cl_idx(end);
hog = zeros(num_cells, args.num_ori_bins);
for i_cl = 1:num_cells
    cl_mask = cl_idx == i_cl;
    g_theta_idx_i = g_theta_idx(cl_mask);
    g_mag_i = g_mag(cl_mask);
    
    hog(i_cl,:) = full(sparse(1,g_theta_idx_i,g_mag_i,1,args.num_ori_bins));
end

switch args.norm_method
    
    case 'l1-sqrt'
        norm = sum(hog(:)) + eps;
        hog = sqrt(hog / norm);
        
    case 'l1-norm'
        norm = sum(hog(:)) + eps;
        hog = hog / norm;
        
    case 'l2-norm'
        norm = sqrt(sum(hog(:).^2) + eps);
        hog = hog / norm;
        
end
        
    
    
    




