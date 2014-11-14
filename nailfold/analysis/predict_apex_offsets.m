function [apex_offset_map apex_class_pred apex_offset_x_pred apex_offset_y_pred] = ...
    predict_apex_offsets(varargin)
%EXTRACT_VESSEL_CENTRES *Insert a one line summary here*
%   [] = extract_vessel_centres()
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
% Created: 18-Jun-2013
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
args = u_packargs(varargin,... % the user's input
    0, ... % non-strict mode
    {'apex_class_rf',...
    'apex_offset_x_rf',...
    'apex_offset_y_rf',...
    'vessel_feature_im',...
    'vessel_centre'},...
    'smoothing_sigma', 2,...
    'num_cells', 8,...
    'hog_args', [],...
    'xy', [],...
    'apex_class_thresh', 0.5,...
    'separate_trees', 0,...
    'max_size', 1000,...
    'base_width', 20,...
    'include_pts', []);
clear varargin;

if ~isempty(args.xy)
    patch_sz = sqrt(size(args.xy,1));
else
    patch_sz = args.num_cells*args.hog_args.cell_sz(1);
    patch_sz = patch_sz + 2; %Account for padding
    patch_sz2 = (patch_sz - 1)/2;

    %Set up x,y coordinates for template patch
    x = repmat(-patch_sz2:patch_sz2, patch_sz, 1);
    y = x';
    args.xy = [x(:) y(:)];
end
hog_sz = args.hog_args.num_ori_bins * args.num_cells^2;
num_pts = length(args.vessel_centre.x);

if args.separate_trees
    num_preds = length(args.apex_offset_x_rf.trees);
else
    num_preds = 1;
end

if isempty(args.include_pts)
    include_pts = true(num_pts,1);
else
    include_pts = args.include_pts;
end
apex_class_pred = zeros(num_pts,1);
apex_offset_x_pred = zeros(num_pts,num_preds);
apex_offset_y_pred = zeros(num_pts,num_preds);

%Create smoothing kernel for feature image
if args.smoothing_sigma 
    if length(args.smoothing_sigma) == 1
        g = gaussian_filters_1d(args.smoothing_sigma);
        g = g / sum(g);
    else
        g = args.smoothing_sigma;
    end
    
    %Smooth the vessel probs
    args.vessel_feature_im = conv2(g', g, args.vessel_feature_im, 'same');
end

curr_pt = 1;
for i_pt = 1:num_pts

    %If the first point in a new part, make a container for the HoGs
    if curr_pt == 1
        part_sz = min(args.max_size, num_pts-i_pt+1);
        vessel_hog = zeros(part_sz, hog_sz);
        start_i = i_pt;
    end

    %Get predicted scale and orientation at this point
    vxc = args.vessel_centre.x(i_pt);
    vyc = args.vessel_centre.y(i_pt);
    ori_c = angle(args.vessel_centre.ori(i_pt))/2;
    width_c = args.vessel_centre.width(i_pt);

    %Get scale relative to base width a make rotation matrix
    rot = [cos(ori_c) -sin(ori_c); sin(ori_c) cos(ori_c)];
    scale = width_c / args.base_width;

    %Transform points given scale and angle and translate to
    %candidate position
    xya = args.xy * rot * scale;
    xa = reshape(xya(:,1) + vxc, [patch_sz patch_sz]);
    ya = reshape(xya(:,2) + vyc, [patch_sz patch_sz]);

    %Sample vessel prob patch
    vessel_feature_patch = interp2(args.vessel_feature_im, xa, ya, '*linear', 0);
    [hog] = compute_HoG(vessel_feature_patch, args.hog_args);       
    vessel_hog(curr_pt,:) = hog(:)';

    %Check if we've filled the part. If so, set the curr_pt back to 1, save 
    % the HoGs, and increment the part number. Other increment curr_pt
    if curr_pt == part_sz
        curr_pt = 1;
        part_idx = start_i:i_pt;

        %Classify whether points in part can point to an apex
        [~,votes] = random_forest_class_predict(args.apex_class_rf, vessel_hog);
        apex_class_pred_i = votes(:,2) / length(args.apex_class_rf.trees); clear votes;

        %Threshold these class scores
        include_pts_i = apex_class_pred_i > args.apex_class_thresh;

        %Save these predictions in the main containers
        apex_class_pred(part_idx,:) = apex_class_pred_i; 
        include_pts(part_idx,:) = include_pts(part_idx,:) & include_pts_i; 

        %Now predict the apex offsets for points above the threshold
        if args.separate_trees
            [~, apex_offset_x_pred(part_idx(include_pts_i),:)] = ...
                random_forest_reg_predict(args.apex_offset_x_rf, vessel_hog(include_pts_i,:));
            [~, apex_offset_y_pred(part_idx(include_pts_i),:)] = ...
                random_forest_reg_predict(args.apex_offset_y_rf, vessel_hog(include_pts_i,:));
        else
            apex_offset_x_pred(part_idx(include_pts_i),:) = ...
                random_forest_reg_predict(args.apex_offset_x_rf, vessel_hog(include_pts_i,:));
            apex_offset_y_pred(part_idx(include_pts_i),:) = ...
                random_forest_reg_predict(args.apex_offset_y_rf, vessel_hog(include_pts_i,:));
        end
    else
        curr_pt = curr_pt + 1;
    end
end
[nrows, ncols] = size(args.vessel_feature_im);

args = rmfield(args, {'vessel_feature_im', 'apex_class_rf', 'apex_offset_x_rf', 'apex_offset_y_rf'});

%Now we can loop through all these points again and update the apex
%location map using the offset predictions and weightings
apex_offset_map = zeros(nrows, ncols);
for i_pt = 1:num_pts
    if include_pts(i_pt)
        vxc = args.vessel_centre.x(i_pt);
        vyc = args.vessel_centre.y(i_pt);
        ori_c = angle(args.vessel_centre.ori(i_pt))/2;
        scale = args.vessel_centre.width(i_pt) / args.base_width;

        %Get scale relative to base width a make rotation matrix
        rot = [cos(ori_c) -sin(ori_c); sin(ori_c) cos(ori_c)];

        %Translate predicted offset from the frame coordinates into real
        %image coordinates
        for i_pred = 1:num_preds
            apex_xy_t = [apex_offset_x_pred(i_pt,i_pred) apex_offset_y_pred(i_pt,i_pred)];
            apex_xy_p = scale*apex_xy_t*rot + [vxc vyc];                

            %Add a blurred vote to the image (width of blur depends on scale of
            %frame)
            [im_r, im_c, gxy] = add_gaussian_blur(apex_xy_p(1), apex_xy_p(2), 2*scale, nrows, ncols);
            if ~isempty(gxy)
                apex_offset_map(im_r, im_c) = apex_offset_map(im_r, im_c) + apex_class_pred(i_pt,1)*gxy;
            end
        end
    end
end          
    
    
