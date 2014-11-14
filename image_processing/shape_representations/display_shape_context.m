function [sc_image_c sc_image sc_outer_mask] = display_shape_context(sc, min_r, output_type)
%display_shape_context *Insert a one line summary here*
%   [masks] = display_shape_context(num_r, min_r, num_theta)
%
% Inputs:
%      num_r - *Insert description of input variable here*
%
%      min_r - *Insert description of input variable here*
%
%      num_theta - *Insert description of input variable here*
%
%
% Outputs:
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 07-Jun-2013
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester
    
[num_r, num_theta] = size(sc);

%Workout sizes of each scale bin
[sc_masks] = make_shape_context_mask(num_r, min_r, num_theta, 0, 1.5);

%Pre-allocate sc_image and outer mask
sc_sz = size(sc_masks,1);
sc_image = zeros(sc_sz, sc_sz);
sc_outer_mask = true(sc_sz, sc_sz);

for i_r = 1:num_r
    for i_theta = 1:num_theta
        
        sc_image(sc_masks(:,:,i_r,i_theta)) = sc(i_r,i_theta);
        sc_outer_mask(sc_masks(:,:,i_r,i_theta)) = 0;
    end
end

if ~exist('output_type', 'var'); output_type = 'c'; end

switch output_type
    
    case 'c'
    
        sc_outer_mask3 = cat(3, sc_outer_mask, sc_outer_mask, sc_outer_mask);
        sc_image_c = complex2rgb(sc_image);
        sc_image_c(sc_outer_mask3) = 1;
        
    case 'p'
        
        sc_outer_mask3 = cat(3, sc_outer_mask, sc_outer_mask, sc_outer_mask);
        z_mask = ~abs(sc_image);
        sc_image_c = sc_image ./ abs(sc_image);
        sc_image_c(z_mask) = 0;
        sc_image_c = complex2rgb(sc_image_c);
        sc_image_c(sc_outer_mask3) = 1;
        
    case 'm'
        sc_image_c = abs(sc_image);
        sc_image_c(sc_outer_mask) = inf;
        
end
        
                    


