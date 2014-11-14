function [resampling_map] = make_resampling_map(prediction_image, fg_mask, fov_mask, gt_map, error_ratio, output_type)
%MAKE_SAMPLED_MAPS *Insert a one line summary here*
%   [] = make_sampled_maps(image_list, sampled_data, save_dir)
%
% Inputs:
%      image_list - *Insert description of input variable here*
%
%      sampled_pts_list - *Insert description of input variable here*
%
%      save_dir - *Insert description of input variable here*
%
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 30-Aug-2012
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester

%Check if we need to thin the mask
switch output_type    
    case {'centre_orientation', 'centre_detection', 'width'}
        fg_mask = bwmorph(fg_mask, 'thin', 'inf');
end
if ~isempty(fov_mask)
    fg_mask(~fov_mask) = 0;
end

%Compute errors for the various output types
switch output_type
    
    case {'detection', 'centre_detection'}
        fg_error = 1 - prediction_image;
        bg_error = prediction_image;
        
    case {'orientation','centre_orientation'}
        %Compute absolute error in prediction and scale between 0 and 1
        fg_error = 2*abs(ori_error(gt_map, prediction_image)) / pi;
        
    case 'width'
        fg_error = 0;
        
end

%Weight the resampling maps between error based sampling an uniform
%sampling
resampling_map = 1 - error_ratio + error_ratio*fg_error;
resampling_map(~fg_mask) = 0;

%For detection also compute the weight background resampling map
switch output_type
    case {'detection', 'centre_detection'}
        bg_resampling_map = 1 - error_ratio + error_ratio*bg_error;
        bg_resampling_map(fg_mask & ~fov_mask) = 0;        
        resampling_map(:,:,2) = bg_resampling_map;
end
        
        


        


    


