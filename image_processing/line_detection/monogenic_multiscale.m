function [energy phase orientation scale] = monogenic_multiscale(image_in, num_scales, min_wave_length, mult, sigma_onf)
%MONOGENIC_PHASE_CONG *Insert a one line summary here*
%   [image_out,phase_cong] = monogenic_phase_cong(image_in,num_scales,min_wave_length,mult,sigma_onf,phase_lims)
%
% Inputs:
%      image_in- *Insert description of input variable here*
%
%      num_scales- *Insert description of input variable here*
%
%      min_wave_length- *Insert description of input variable here*
%
%      mult- *Insert description of input variable here*
%
%      sigma_onf- *Insert description of input variable here*
%
%      feature_type- *Insert description of input variable here*
%
%
% Outputs:
%      image_out- *Insert description of input variable here*
%
%      phase_cong- *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 21-Apr-2010
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester 

%Compute monogenic signal of image
[local_amp, local_phase, local_ori] = monogenic(image_in, num_scales, min_wave_length, mult, sigma_onf);

%Get energy as maximum of image amplitudes and the scale associated with
%each maxima
[energy scale] = max(local_amp(:,:,2:num_scales+1), [], 3);

%Loop through levels, taking the phase and orientations from pixels with maximal
%amplitudes at that scale
[r c] = size(image_in);
phase = zeros(r, c);
orientation = zeros(r, c);
for level = 1:num_scales
    scale_idx = scale == level;
    
    ori_scale = local_ori(:,:,level+1);
    orientation(scale_idx) = ori_scale(scale_idx);
    
    phase_scale = local_phase(:,:,level+1);
    phase(scale_idx) = phase_scale(scale_idx);

end

orientation = mod(orientation, pi);

