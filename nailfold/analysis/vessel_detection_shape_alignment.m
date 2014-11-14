function [unaligned_shape_data aligned_shape_data] = vessel_detection_shape_alignment(varargin)
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
args = u_packargs(varargin, '0', {... 
    'contour_files',...
    'vessel_files' },...
    'model_dir', [],...
    'num_shape_pts', 31,...
    'plot', 0);
clear varargin;

num_vessels = length(args.contour_files);

%--------------------------------------------------------------------------
% 1) Extract the original shapes from the data files
apex_shapes = zeros(num_vessels, 2*args.num_shape_pts);
apex_areas = zeros(num_vessels,1);
apex_widths = zeros(num_vessels,1);

for i_v = 1:num_vessels
    %load data
    contour_struc = load(args.contour_files{i_v});
    vessel_struc = u_load(args.vessel_files{i_v});
    vessel_top = vessel_struc.v_pts_new;

    %Save shape into main container and compute area enclosed by shape
    apex_shapes(i_v,:) = vessel_top(:);
    apex_areas(i_v) = polyarea(vessel_top(:,1), vessel_top(:,2));

    %Work out the width of the vessel at the point closest to the vessel
    %apex
    vessel_centre = (contour_struc.outer_edge + contour_struc.inner_edge)/2;
    v_widths_i = sqrt(sum((contour_struc.inner_edge - contour_struc.outer_edge).^2,2));
    dists = sum(bsxfun(@minus, vessel_centre, vessel_top(16,:)).^2,2);
    [~, apex_i] = min(dists);
    apex_widths(i_v) = v_widths_i(apex_i);

end
unaligned_shape_data.apex_shapes = apex_shapes;
unaligned_shape_data.apex_widths = apex_widths;
unaligned_shape_data.apex_areas = apex_areas;
unaligned_shape_data.apex_names = args.vessel_files;

save([args.model_dir 'unaligned_apex_shape_data.mat'], 'unaligned_shape_data');

%--------------------------------------------------------------------------
% 2) Align the shapes to create a mean target apex

%Align shapes
[a_shapes, a_scales, ~, a_rots, a_trans] = ...
    align_shapes(apex_shapes, 'area', mean(apex_areas), 'shiftOrigin', 0);

%Get mean of aligned shape
mean_shape = [mean(a_shapes(:,1:args.num_shape_pts))',...
                mean(a_shapes(:,args.num_shape_pts+1:2*args.num_shape_pts))'];

%Compute the rotation needed so that the mean stands vertically (i.e. the vector from
%its centre of mass to the apex is at 90 degrees)
mean_centre = mean(mean_shape);
mean_vec = mean_centre - mean_shape(16,:);
mean_theta = atan2(mean_vec(2), mean_vec(1));
apex_theta = rotation_matrix(3*pi/2+mean_theta);

%Now rotate each aligned shape so that it stands vertically - for shapes
%that were reflected in their initial alignment, undo the reflection
for i_sh = 1:size(a_shapes,1)
    a_shape_i = reshape(a_shapes(i_sh,:), [], 2);
    a_shape_i = a_shape_i * apex_theta;
    if sign(a_rots(1,1,i_sh)) ~= sign(a_rots(2,2,i_sh))        
        a_shape_i(:,1) = -a_shape_i(:,1);
        
        if args.plot
            figure; axis ij equal; hold all;
            plot(-a_shape_i(:,1), a_shape_i(:,2));
            plot(a_shape_i(:,1), a_shape_i(:,2));
        end
    end
    
    %Save the realigned shape
    a_shapes(i_sh,:) = a_shape_i(:)';

    %Now compute the transform needed to map each original shape to its
    %aligned counterpart
    u_shape_i = reshape(apex_shapes(i_sh,:), [], 2);
    [dd, ~, t] = mb_procrustes(a_shape_i, u_shape_i, 0);
    a_scales(i_sh) = t.b;
    a_rots(:,:,i_sh) = t.T;
    a_trans(i_sh,:) = t.c(1,:) / t.b;
    
    %Note this alignment should be perfect, so check dd is neglible
    if dd > 1e-5
        error('Shape did not align correctly');
    end
end

%Compute the new mean shape
aligned_shape_data.a_shapes = a_shapes;
aligned_shape_data.a_scales = a_scales;
aligned_shape_data.a_rots = a_rots;
aligned_shape_data.a_trans = a_trans;
save([args.model_dir 'aligned_apex_shape_data.mat'], 'aligned_shape_data');

if args.plot
    mean_shape = reshape(mean(a_shapes), [], 2);
    figure; axis ij equal; hold all;
    plot(a_shapes(:,1:31)', a_shapes(:,32:62)'); 
    plot(mean_shape(:,1), mean_shape(:,2), 'k', 'linewidth', 2)
end
