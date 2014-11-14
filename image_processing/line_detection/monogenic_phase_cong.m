function [image_out phase_cong orientation local_ori] = monogenic_phase_cong(image_in, num_scales, min_wave_length, mult, sigma_onf, feature_type)
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
if nargin < 6
    feature_type = 'pos_line';
end

%Compute monogenic signal of image
[local_amp, local_phase, local_ori] = monogenic(image_in, num_scales, min_wave_length, mult, sigma_onf);

%Loop through levels summing local phase/amp complex numbers at ecah pixel
phase_cong = 0;
for level = 2:num_scales+1
    phase_cong = phase_cong + local_amp(:,:,level).*exp(i*local_phase(:,:,level));
end

%Now return the magnitude
switch feature_type
    case 'pos_line'
        image_out = imag(phase_cong);
        image_out(image_out < 0) = 0;
        
    case 'neg_line'
        image_out = imag(phase_cong);
        image_out(image_out > 0) = 0;
        
    case 'edge'
        image_out = real(phase_cong);
    otherwise
        %Warn that feature type not recoginsed... process as line
        image_out = real(phase_cong);
end
image_out = image_out / max(image_out(:));

if nargout > 2
    [dummy max_scale_idx] = max(local_amp(:,:,2:end),[],3);
    orientation = zeros(size(image_in));
    for level = 2:num_scales+1
        level_oris = local_ori(:,:,level);
        orientation(max_scale_idx==level-1) = level_oris(max_scale_idx==level-1); 
    end
    orientation = 180*orientation / pi;
    orientation(orientation < 0) = orientation(orientation < 0) + 180;
end