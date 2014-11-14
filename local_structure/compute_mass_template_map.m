function [template_map scales_map] = compute_mass_template_map(mammo, radii_scales, radii_ratio)
%COMPUTE_MASS_TEMPLATE_MAP *Insert a one line summary here*
%   [] = compute_mass_template_map()
%
% Inputs:
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 22-Oct-2010
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

%pre-allocate template and scale maps
template_map = -ones(size(mammo));
scales_map = uint8(zeros(size(mammo)));

%Go through each scale template matching
for ii = 1:length(radii_scales)
    
    %Create template
    R = radii_scales(ii);
    R2 = round(radii_ratio*R);

    xx = repmat(1:2*R2, 2*R2, 1) - R2;
    yy = xx';
    template = R^2 - xx.^2 - yy.^2;
    template(xx.^2 + yy.^2 > R^2) = 0;

    %Do template matching
    temp_template = normxcorr2(template, double(mammo));
    temp_template = temp_template(R2:end-R2,R2:end-R2);

    %Work out which pixels are now maximal for this scale
    scale_mask = temp_template > template_map;
    
    %Record template matching score and scale for those pixels
    template_map(scale_mask) = temp_template(scale_mask);
    scales_map(scale_mask) = R;
end