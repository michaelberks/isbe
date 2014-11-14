function [] = compute_mosiac_midline_set(fov_mask_dir, output_dir, varargin)
%COMPUTE_MOSIAC_MIDLINE *Insert a one line summary here*
%   [centres_r, rot_mat, ncolsr, nrowsr, ncols, nrows, bad] = compute_mosiac_midline(fov_mask, poly_n)
%
% Inputs:
%      fov_mask - *Insert description of input variable here*
%
%      poly_n - *Insert description of input variable here*
%
%
% Outputs:
%      centres_r - *Insert description of input variable here*
%
%      rot_mat - *Insert description of input variable here*
%
%      ncolsr - *Insert description of input variable here*
%
%      nrowsr - *Insert description of input variable here*
%
%      ncols - *Insert description of input variable here*
%
%      nrows - *Insert description of input variable here*
%
%      bad - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 04-Dec-2013
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

args = u_packargs(varargin,... % the user's input
    0, ... % non-strict mode
    'start_i', 1,...
    'end_i', [],...
    'ori_bins', -45:45,...
    'resize_factor', 8, ...
    'bad_mosaic_sz', 55,...
    'border_size', 2,...
    'min_depth', 470,...
    'poly_n', 5);
clear varargin;

im_list = dir([fov_mask_dir '*.mat']);
if isempty(args.end_i)
    args.end_i = length(im_list);
end
create_folder(output_dir);

display(['Processing images ' num2str(args.start_i) ' to ' num2str(args.end_i)]);
for i_im = args.start_i:args.end_i
    display(['Processing image : ' num2str(i_im)]);
    
    im_name = im_list(i_im).name(1:6);
    fov_mask = u_load([fov_mask_dir im_name '_f_mask.mat']);

    [centres_r, rot_mat, ncolsr, nrowsr, ncols, nrows, bad] =...
        compute_mosiac_midline(fov_mask, args); %#ok
    
    save([output_dir im_name '.mat'], 'centres_r', 'rot_mat', 'ncolsr', 'nrowsr', 'ncols', 'nrows', 'bad');
    
end


