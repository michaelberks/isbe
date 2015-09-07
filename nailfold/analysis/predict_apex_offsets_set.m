function [] = predict_apex_offsets_set(varargin)
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
    'task_id',              unixenv('SGE_TASK_ID',1), ...
    'num_jobs',             unixenv('NUM_JOBS',100), ...
    'model_id',             unixenv('MODEL_PATH', []),...
    'model_root',           [unixenv('DATA_ROOT',[]) unixenv('MODEL_ROOT',[nailfoldroot,'models/apex'])], ...
    'model_name',           unixenv('MODEL_NAME', 'predictor'),...
    'data_dir',             [unixenv('DATA_ROOT',[]) unixenv('IMAGE_ROOT', [nailfoldroot 'data/rsa_study/test'])],...
    'feature_im_dir',       unixenv('IMAGE_DIR', 'predictions/detection/rf_classification/257273'),...
    'feature_im_ext',       '*.mat',...
    'fov_mask_dir',         unixenv('MASK_DIR',[]),...
    'centre_dir',           unixenv('FG_MASK_DIR', 'vessel_centres'),...
    'apex_map_dir',         'apex_maps',...
    'max_size',             unixenv('MAX_SIZE', 1000),...
    'separate_trees',       0,...
    'smoothing_sigma',      1,...
    'num_cells',            8,...
    'cell_sz',              8,... %Size of HoG cells in blocks
    'block_sz',             [2 2],...%Size of blocks in cells
    'num_ori_bins',         9,... %Number of bins in orientation histograms
    'norm_method',          'none',... %Method for local normalisation
    'block_spacing',        8,... %Separation between blocks in pixels - controls level of overlap
    'gradient_operator',    [-1 0 1],...
    'spatial_sigma',        0, ...
    'angle_wrap',           1,...
    'base_width',           20,...
    'apex_class_thresh',    unixenv('THRESH', 0),...
    'overwrite',            unixenv('OVERWRITE',false));
clear varargin;

%Load in the classification and offset models
rf_class_dir = [args.model_root '/classification/' args.model_id '/'];
rf_offset_x_dir = [args.model_root '/offset_x/' args.model_id '/'];
rf_offset_y_dir = [args.model_root '/offset_y/' args.model_id '/'];
apex_class_rf = u_load([rf_class_dir args.model_name '.mat']);
apex_offset_x_rf = u_load([rf_offset_x_dir args.model_name '.mat']);
apex_offset_y_rf = u_load([rf_offset_y_dir args.model_name '.mat']);

%Form full directory paths and create folder for HoGs
fov_mask_dir = [args.data_dir '/' args.fov_mask_dir '/'];
feature_im_dir = [args.data_dir '/' args.feature_im_dir '/'];
centre_dir = [args.data_dir '/' args.centre_dir '/'];
apex_map_dir = [args.data_dir '/' args.apex_map_dir '/' args.model_id '/'];
create_folder(apex_map_dir);

%Get HoG args from main args
hog_args.cell_sz = [args.cell_sz args.cell_sz];
hog_args.block_sz = args.block_sz;
hog_args.num_ori_bins = args.num_ori_bins;
hog_args.norm_method = args.norm_method;
hog_args.block_spacing = args.block_spacing;
hog_args.gradient_operator = args.gradient_operator;
hog_args.spatial_sigma = args.spatial_sigma;
hog_args.angle_wrap = args.angle_wrap;

%1. Workout number of images in job
image_list = dir([feature_im_dir args.feature_im_ext]);
num_images = length(image_list);
job_size = ceil(num_images / args.num_jobs);

%2. Workout start and end indices for job
start_i	= (args.task_id-1)*job_size + 1;
end_i	= min(args.task_id*job_size, num_images);
display(['Computing HoGs for images ' num2str(start_i) ' to ' num2str(end_i)]);

%Create smoothing kernel for feature image
if args.smoothing_sigma
    g = gaussian_filters_1d(args.smoothing_sigma);
    g = g / sum(g);
end

%Get patch size and form template x,y coordinates for the patch
patch_sz = args.num_cells*args.cell_sz;
patch_sz = patch_sz + 2; %Account for padding
patch_sz2 = (patch_sz - 1)/2;

%Get hog size from the hog_args struct
%Set up x,y coordinates for template patch
x = repmat(-patch_sz2:patch_sz2, patch_sz, 1);
y = x';
xy = [x(:) y(:)];

%Loop though each image
for i_im = start_i:end_i
    im_name = image_list(i_im).name(1:end-9);
    
    save_name = [apex_map_dir image_list(i_im).name(1:end-3) 'mat'];
    
    if ~args.overwrite && exist(save_name, 'file')
        display(['Image ' num2str(i_im) ...
            ' already exists. Skipping (switch overwrite to 1 if necessary)']);
        continue;
    end
    
    display(['Processing image ' num2str(i_im) ', ' datestr(now)]);  
    
    %Load in images and vessel markup
    vessel_feature_im = u_load([feature_im_dir image_list(i_im).name]);
    load([centre_dir im_name '_vc.mat'], 'vessel_centre', 'nrows', 'ncols');
    
    %If we've been given a mask discard pts from the edges of the mask
    if ~isempty(args.fov_mask_dir)
        display('Using mask');
        f_mask = u_load([fov_mask_dir im_name '_f_mask.mat']);
        [discard_pts] = discard_edge_preds(vessel_centre, f_mask);
        include_pts = ~discard_pts;
    else
        include_pts = [];
    end

    [apex_offset_map apex_class_pred apex_offset_x_pred apex_offset_y_pred] = ...
        predict_apex_offsets(...
            'apex_class_rf', apex_class_rf,...
            'apex_offset_x_rf', apex_offset_x_rf,...
            'apex_offset_y_rf', apex_offset_y_rf,...
            'vessel_feature_im', vessel_feature_im, ...
            'vessel_centre', vessel_centre, ...
            'separate_trees', args.separate_trees,...
            'smoothing_sigma', g,...
            'num_cells', args.num_cells,...
            'hog_args', hog_args,...
            'xy', xy,...
            'apex_class_thresh', args.apex_class_thresh,...
            'max_size', args.max_size,...
            'base_width', args.base_width,...
            'include_pts', include_pts); %#ok
        
    %Save this map
    save(save_name, 'apex_offset_map',...
        'apex_class_pred', 'apex_offset_x_pred', 'apex_offset_y_pred');        
end          
    
    
