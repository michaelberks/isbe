function [] = compute_apex_candidate_hogs(varargin)
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
args = u_packargs(varargin, 0, ... % the user's input
    {'image_names'},           ...
    'data_dir',             [nailfoldroot 'data/rsa_study/master_set/'],...
    'feature_im_dir',       'images',...
    'ori_dir',              'rf_regression/296621/',...
    'width_dir',            'rf_regression/297037/',...
    'candidates_dir',       'apex_maps/local_maxima',...
    'hog_dir',              'apex_hogs',...
    'feature_sigma',        0,...
    'prob_sigma',           2,...
    'ori_sigma',            0,...
    'width_sigma',          2,...
    'num_cells',            8,...
    'cell_sz',              8,... %Size of HoG cells in blocks
    'block_sz',             [2 2],...%Size of blocks in cells
    'num_ori_bins',         12,... %Number of bins in orientation histograms
    'norm_method',          'l1-sqrt',... %Method for local normalisation
    'block_spacing',        8,... %Separation between blocks in pixels - controls level of overlap
    'gradient_operator',    [-1 0 1],...
    'spatial_sigma',        0, ...
    'angle_wrap',           1,...
    'base_width',           20, ...
    'dist_thresh',          24^2,...
    'overwrite',            0);

%Form full directory paths and create folder for HoGs
feature_im_dir = [args.data_dir args.feature_im_dir '/'];
ori_dir = [args.data_dir 'predictions/orientation/' args.ori_dir '/'];
width_dir = [args.data_dir 'predictions/width/' args.width_dir '/'];
candidates_dir = [args.data_dir '/' args.candidates_dir '/'];
hog_dir = [args.data_dir '/' args.hog_dir '/'];
create_folder(hog_dir);

%Get HoG args from main args
hog_args.cell_sz = [args.cell_sz args.cell_sz];
hog_args.block_sz = args.block_sz;
hog_args.num_ori_bins = args.num_ori_bins;
hog_args.norm_method = args.norm_method;
hog_args.block_spacing = args.block_spacing;
hog_args.gradient_operator = args.gradient_operator;
hog_args.spatial_sigma = args.spatial_sigma;
hog_args.angle_wrap = args.angle_wrap;

num_images = length(args.image_names);


%Get patch size and form template x,y coordinates for the patch
patch_sz = args.num_cells*args.cell_sz;
patch_sz = patch_sz + 2; %Account for padding
patch_sz2 = (patch_sz - 1)/2;
hog_sz = args.num_cells*args.num_cells*args.num_ori_bins;

%Set up x,y coordinates for template patch
x = repmat(-patch_sz2:patch_sz2, patch_sz, 1);
y = x';
xy = [x(:) y(:)];

if args.feature_sigma
    g_feat = gaussian_filters_1d(args.feature_sigma);
    g_feat = g_feat / sum(g_feat);    
end
if args.ori_sigma
    g_ori = gaussian_filters_1d(args.ori_sigma);
    g_ori = g_ori / sum(g_ori);
end 
if args.width_sigma
    g_width = gaussian_filters_1d(args.width_sigma);
    g_width = g_width / sum(g_width);
end 

%Loop though each image
for i_im = 1:num_images
    
    im_name = args.image_names{i_im} ;  
    
    if ~args.overwrite && exist([hog_dir im_name '_hog.mat'], 'file')
        display(['Skipping image ' num2str(i_im) ',' hog_dir im_name '_hog.mat already exists. Re-run with overwrite=1 if necessary']); 
        continue;
    end
    display(['Processing image ' num2str(i_im) ', ' datestr(now)]);  
    
    %Load in images and vessel markup
    vessel_feature_im = double(u_load([feature_im_dir im_name '.mat']));
    vessel_ori = u_load([ori_dir im_name '_pred.mat']);
    vessel_width = u_load([width_dir im_name '_pred.mat']);
    load([candidates_dir im_name '_candidates'], 'candidate_xy');
        
    if isempty(candidate_xy) %#ok
        candidates_hogs = []; %#ok
        candidate_oris = []; %#ok
        candidate_widths = []; %#ok
    else
    
        if args.feature_sigma
            vessel_feature_im = conv2(g_feat', g_feat, vessel_feature_im, 'same');
        end
        if args.ori_sigma
            vessel_ori = conv2(g_ori', g_ori, vessel_ori, 'same');
        end 
        if args.width_sigma
            vessel_width = conv2(g_width', g_width, vessel_width, 'same');
        end
        candidate_oris = interp2(vessel_ori, candidate_xy(:,1), candidate_xy(:,2));
        candidate_widths = interp2(vessel_width, candidate_xy(:,1), candidate_xy(:,2));
        clear vessel_ori vessel_width;

        candidate_oris = angle(candidate_oris / 2); 

        num_candidates = size(candidate_xy,1);
        candidates_hogs = zeros(num_candidates, hog_sz);

        for i_can = 1:num_candidates

            %EXtract patch and compute HoG
            ori_c = candidate_oris(i_can);
            width_c = candidate_widths(i_can);

            %Get scale relative to base width a make rotation matrix
            rot = [cos(ori_c) -sin(ori_c); sin(ori_c) cos(ori_c)];
            scale = width_c / args.base_width;

            %Transform points given scale and angle and translate to
            %candidate position
            xya = xy * rot * scale;
            xa = reshape(xya(:,1) + candidate_xy(i_can,1), patch_sz, patch_sz);
            ya = reshape(xya(:,2) + candidate_xy(i_can,2), patch_sz, patch_sz);

            %Sample vessel prob patch
            vessel_feature_patch = interp2(vessel_feature_im, xa, ya, '*linear', 0);
            [hog] = compute_HoG(vessel_feature_patch, hog_args);       
            candidates_hogs(i_can,:) = hog(:)';
        end
    end
    save([hog_dir im_name '_hog.mat'], 'candidates_hogs', 'candidate_oris', 'candidate_widths');
end          
    
    
