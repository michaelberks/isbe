function [template_scores] = ...
    targetted_template_matching(feature_map, potential_map, scale_map, orientation_map, template, varargin)
%TARGETTED_TEMPLATE_MATCHING *Insert a one line summary here*
%   [maxima_pos, maxima_vals] = targetted_template_matching(varargin)
%
% TARGETTED_TEMPLATE_MATCHING uses the U_PACKARGS interface function
% and so arguments to the function may be passed as name-value pairs
% or as a struct with fields with corresponding names. The field names
% are defined as:
%
% Mandatory Arguments:
%
% Optional Arguments:
%
% Outputs:
%      maxima_pos - *Insert description of input variable here*
%
%      maxima_vals - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 30-Apr-2013
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

% Unpack the arguments:
%args = u_packargs(varargin, '0', {}, {});
%clear varargin;
args = u_packargs(varargin, '0', ...
    'filter_scales_map', [],...
    'plot', 0);
clear varargin;

[rows, cols, ~] = size(feature_map);

%Get indices and subscripts of potential points
potential_idx = find(potential_map);
[potential_r potential_c] = ind2sub([rows cols], potential_idx);

%Sample scales and orientations at the potential points
potential_scales = scale_map(potential_idx);
potential_thetas = orientation_map(potential_idx);

if ~isempty(args.filter_scales_map)
    potential_filter_scales = args.filter_scales_map(potential_idx);
end

%Get template size and set up sampling pts matrix
template_sz = size(template,1);
lim = (template_sz-1)/2;
x = repmat(-lim:lim, template_sz, 1);
y = x';
xy = [x(:) y(:)];

%Pre-allocate output
template_scores = zeros(rows, cols);
size_T = template_sz^2;

template_mu = mean(template(:));
template_sd = std(template(:), 1);
template = template(:)' / size_T;


for i_pt = 1:length(potential_idx)
    
    %Transform the smapling points by the scale and orientation at this
    %location, then shift to the location
    R_theta = rotation_matrix(potential_thetas(i_pt));
    xya = xy * R_theta * potential_scales(i_pt);
    xa = xya(:,1) + potential_c(i_pt);
    ya = xya(:,2) + potential_r(i_pt);
    
    %Sample from the feature map at these point
    if ~isempty(args.filter_scales_map)
        f_scale = potential_filter_scales(i_pt);
        feature_patch = interp2(feature_map(:,:,f_scale), xa, ya, '*bilinear');
    else
        feature_patch = interp2(feature_map, xa, ya, '*bilinear');
    end
    
    feature_patch(isnan(feature_patch)) = 0;
    feature_sum2 = sum(feature_patch.^2);
    feature_mu = sum(feature_patch) / size_T;    
    feature_sd = sqrt( max((feature_sum2/size_T - feature_mu^2),0) ); 
        
    %Compute NCC against template
    template_scores(potential_idx(i_pt)) =...
        ( template*feature_patch - feature_mu*template_mu ) / ...
        (template_sd * feature_sd);
end
    
        
        
    

