function [apex_template] = vessel_detection_template_building(varargin)
%VESSEL_DETECTION_TEMPLATE_BUILDING *Insert a one line summary here*
%   [] = vessel_detection_template_building(varargin)
%
% VESSEL_DETECTION_TEMPLATE_BUILDING uses the U_PACKARGS interface function
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
% Created: 01-May-2013
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

% Unpack the arguments:
args = u_packargs(varargin, '0', { ...
    'apex_shapes',...
    'apex_widths',...
    'a_scales',...
    'a_rots'},...
    'model_dir', [],...
    'vessel_files', [], ...
    'vessel_prob_files', [],...
    'template_size', 49,...
    'do_intensity', 0,...
    'do_g1', 0,...
    'do_g2', 1,...
    'do_vessel_prob', 0,...
    'sigma2', 4,...
    'sigma1', 4,...
    'plot', 0);
clear varargin;

%--------------------------------------------------------------------------
% Now loop through each vessel using the shape transforms to extract a
% scaled and rotated patch about each apex
load_vessel = false;
if args.do_intensity
    apex_template.intensity = zeros(args.template_size, args.template_size);
    load_vessel = true;
end
if args.do_g1
    apex_template.g1 = zeros(args.template_size, args.template_size);
    load_vessel = true;
end
if args.do_g2
    apex_template.g2 = zeros(args.template_size, args.template_size);
    load_vessel = true;
end
if args.do_vessel_prob
    apex_template.vessel_prob = zeros(args.template_size, args.template_size);
end    
    
%Make initial sampling points
lim = (args.template_size-1)/2;
x = repmat(-lim:lim, 49, 1);
y = x';
xy = [x(:) y(:)];

%Compute mean scale
mean_scale = mean(args.a_scales);

[num_vessels num_shape_pts] = size(args.apex_shapes);
num_shape_pts = num_shape_pts/2;
top_ix = (num_shape_pts + 1) / 2;
top_iy = num_shape_pts + top_ix;

for i_v = 1:num_vessels
    
    %Scale and rotate the sampling points
    scale_factor = mean_scale / args.a_scales(i_v);        
    xya = (xy * args.a_rots(:,:,i_v)' * scale_factor);
    xa = reshape(xya(:,1) + args.apex_shapes(i_v,top_ix), args.template_size, args.template_size);
    ya = reshape(xya(:,2) + args.apex_shapes(i_v,top_iy), args.template_size, args.template_size);
    
    %If needed, load vessel patch
    if load_vessel || (args.plot && i_v <= 5)
        vessel_struc = u_load(args.vessel_files{i_v});
        vessel_patch = equalise_nailfold_intensity(vessel_struc.vessel_patch);
    end

    %Sample directly from vessel patch
    if args.do_intensity
        intensity_patch = interp2(vessel_patch, xa, ya, 'bilinear');
        intensity_patch(isnan(intensity_patch)) = mean(intensity_patch(~isnan(intensity_patch)));
        apex_template.intensity(:,:) = apex_template.intensity(:,:) +...
            intensity_patch;
    end
    
    %Compute Gaussian 1st derivatives and sample from patch
    if args.do_g1
        [mag_1d] = gaussian_1st_derivative_gradient(vessel_patch, args.sigma1*scale_factor);
        apex_template.g1(:,:) = apex_template.g1(:,:) +...
            interp2(mag_1d, xa, ya, 'bilinear', 0);
    end
    
    %Compute Gaussian 2nd derivatives and sample from patch
    if args.do_g2    
        [mag_2d] = gaussian_2nd_derivative_line(vessel_patch, args.sigma2*scale_factor);
        apex_template.g2(:,:) = apex_template.g2(:,:) +...
            interp2(mag_2d, xa, ya, 'bilinear', 0);
    end
    
    %Load vessel probability map and sample
    if args.do_vessel_prob
        vessel_prob = u_load(args.vessel_prob_files{i_v});
        vessel_prob_patch = interp2(vessel_prob, xa, ya, 'bilinear');
        vessel_prob_patch(isnan(vessel_prob_patch)) = 1;
        
        apex_template.vessel_prob(:,:) = apex_template.vessel_prob(:,:) +...
            vessel_prob_patch;
    end

    if args.plot && i_v <= 5
        figure; imgray(vessel_patch);
        plot(xa, ya, 'b.', 'markersize', 3);
        plot(...
            args.apex_shapes(i_v,1:num_shape_pts),...
            args.apex_shapes(i_v,num_shape_pts+1:2*num_shape_pts), 'y');
    end
end
save([args.model_dir 'apex_templates.mat'], 'apex_template');