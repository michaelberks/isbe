function [image_out] = ...
	add_image_noise(image_in, noise_type, noise_params)
%ADD_IMAGE_NOISE Adds noise to an image given some noise model
%   [image_out] = add_image_noise(image_in, noise_type, noise_params)
%
% Inputs:
%      
%      image_in - input image (noise may be signal dependant)
% 
%      noise_type - name of noise {'gaussian', rician'}
% 
%      noise_params - parameters to control the noise
%
% Outputs:
%      image_out - Resulting image
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 02-Feb-2011
% Author: Phil Tresadern 
% Email : philip.tresadern@manchester.ac.uk 
% Phone : +44 (0)161 275 5114 
% Copyright: (C) University of Manchester 

switch noise_type
    
    case 'gaussian'
        [sz_y sz_x] = size(image_in);
        image_out = image_in + reshape(sample_from_normal(0, noise_params(1)^2, sz_x*sz_y), sz_y, sz_x);
        
    case 'rician'
        
        image_out = ricernd(image_in, noise_params(1));
        
    %TO DO: other noise - poisson etc.
        
    otherwise
        error(['Noise type: ' noise_type ' not recognised']);
end