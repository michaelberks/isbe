function mean_orientation_map =...
    build_orientation_model(mean_shape, mammo_names, varargin)
%WARP_MASS_CENTRES *Insert a one line summary here*
%   [mass_centres] = warp_mass_centres(mean_shape,mammo_names,args.plot)
%
% Inputs:
%      mean_shape- *Insert description of input variable here*
%
%      mammo_names- *Insert description of input variable here*
%
%      args.plot- *Insert description of input variable here*
%
%
% Outputs:
%      mass_centres- *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 17-Nov-2009
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester 

% Unpack the arguments:
args = u_packargs(varargin,... % the user's input
    '0', ... % non-strict mode
    'mlo', 0,...
    'seg_dir', 'C:\isbe\asymmetry_project\data\segmentations\2004_screening\normals\',...
    'ori_dir', 'C:\isbe\asymmetry_project\data\orientation_maps\rf_small\2004_screening_processed\normals\',...
    'sigma', 0,...
    'plot', 0);
clear varargin;

%Mean breast mask
mean_mask = poly2mask(...
        mean_shape(:,1) - min(mean_shape(:,1)),...
        mean_shape(:,2) - min(mean_shape(:,2)),...
        ceil(max(mean_shape(:,2)) -  min(mean_shape(:,2))),...
        ceil(max(mean_shape(:,1)) -  min(mean_shape(:,1))));
    
%Get xy points inside mean breast
[mean_y mean_x] = find(mean_mask);
mean_x = mean_x + min(mean_shape(:,1));
mean_y = mean_y + min(mean_shape(:,2));

%Set up mean orientation map
mean_orientation_map = zeros(size(mean_mask));

%Get number of masses
n_mammos = length(mammo_names);

%Warp each breast shape to the breast mean, and in doing so compute the
%warped value of the mass centroid
for ii = 1:n_mammos
    
    %Load in segmentation
    seg_list = dir([args.seg_dir '*' mammo_names{ii} '*.mat']);
    seg = u_load([args.seg_dir seg_list(1).name]);
    
    %Load in orientation_mask
    ori_list = dir([args.ori_dir '*' mammo_names{ii} '*.mat']);
    orientation_map = u_load([args.ori_dir ori_list(1).name]);
    
    %smooth map if necessary
    if args.sigma
        orientation_map = imfilter(orientation_map, ...
            fspecial('gaussian', 5*args.sigma+1, args.sigma));
    end
    %orientation_map = exp(i*2*angle(orientation_map));
    %orientation_map = imresize(orientation_map, seg.size);
    
    %Make mask of breast
    breast_mask = poly2mask(...
        seg.breast_border(:,1),...
        seg.breast_border(:,2),...
        seg.size(1),...
        seg.size(2));
    
    %Get xy points inside breast
    if mammo_names{ii}(4) == 'R'
        breast_mask = fliplr(breast_mask);
        orientation_map = fliplr(orientation_map);
        orientation_map = -orientation_map;    
    end
    [breast_y breast_x] = find(breast_mask);
    breast_xy = [breast_x breast_y]; clear breast_x breast_y;
    
    %Get breast border as a standard CC/MLO shape
    if args.mlo
        [breast_shape dummy T] = get_mlo_shape(seg_list, args.seg_dir, 50);
        breast_shape = reshape(breast_shape, [], 2);
        breast_xy = breast_xy*T.rot - repmat(T.t,size(breast_xy,1),1);

    else
        breast_shape = reshape(get_cc_shape(seg_list, args.seg_dir, 50),[],2);
    end    
    
    %Align breast shape to the target mean shape
    [dd breast_shape_a t] = mb_procrustes(mean_shape, breast_shape);
    
    %Map mass border and mass centre into the space of the aligned shape
    breast_xy_a = t.b*breast_xy*t.T + repmat(t.c(1,:), size(breast_xy,1),1);
    
    %Now thin-plate spline warp each breast shape to the target mean shape
    
    %Breast shape form the src points
    s_x = breast_shape_a(:,1)';
    s_y = breast_shape_a(:,2)';
    
    %Mean shape forms the target displacements
    z_x = mean_shape(:,1)';
    z_y = mean_shape(:,2)';
    
    %Build geometric transform
    T = geom_alignpoints_2d([z_x; z_y], [s_x; s_y], [],...
        'transform', 'spline');    
    
    %Warp the breast shape, mass border and mass centre to the mean shape
    %space
    [breast_xy_w] = geom_transformpoints(breast_xy_a', T)';
    [breast_shape_w] = geom_transformpoints(breast_shape_a', T)';        
    clear T;
    
    %Now sample orientations in the mean mask by interpolating the warped
    %pts
    mean_orientation_map(mean_mask) = mean_orientation_map(mean_mask) + griddata(...
        breast_xy_w(:,1),...
        breast_xy_w(:,2),...
        orientation_map(breast_mask),...
        mean_x, mean_y, 'nearest');
    
    %If debugging display aligned and warped breast shape and mass
    if args.plot
        figure; 
        subplot(1,2,1); hold on;
        plot(mean_shape(:,1), mean_shape(:,2), 'b');
        plot(breast_shape_w(:,1), breast_shape_w(:,2), 'r:');
        plot(breast_xy_w(1:64:end,1), breast_xy_w(1:64:end,2), 'c.');
        axis equal;
        subplot(1,2,2); image(complex2rgb(mean_orientation_map)); axis image;
    end
    clear breast_mask orientation_map breast_xy*;
end