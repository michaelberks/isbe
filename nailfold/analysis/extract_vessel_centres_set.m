function [] = extract_vessel_centres_set(varargin)
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
    '0', ... % non-strict mode
    'task_id',              unixenv('SGE_TASK_ID',1), ...
    'num_jobs',             unixenv('NUM_JOBS',100), ...
    'data_dir',             [nailfoldroot 'data/rsa_study/test'],...
    'prob_dir',             'rf_classification/257273',...
    'ori_dir',              'rf_regression/259076',...
    'width_dir',            'rf_regression/257847',...
    'fov_mask_dir',         'fov_masks',...
    'centre_dir',           'vessel_centres',...
    'prob_sigma',           2,...
    'ori_sigma',            0,...
    'width_sigma',          2,...
    'strong_vessel_thresh', 0.25,...
    'weak_vessel_thresh',   0,...
    'overwrite',            0);

prob_dir = [args.data_dir 'predictions/detection/' args.prob_dir '/'];
ori_dir = [args.data_dir 'predictions/orientation/' args.ori_dir '/'];
width_dir = [args.data_dir 'predictions/width/' args.width_dir '/'];
centre_dir = [args.data_dir '/' args.centre_dir '/'];
create_folder(centre_dir);

if isempty(args.fov_mask_dir)
    fov_mask_dir = [];
else
    fov_mask_dir = [args.data_dir '/' args.fov_mask_dir '/'];
end

pred_list = dir([prob_dir '*.mat']);
num_images = length(pred_list);

%1. Workout number of images in job
job_size = ceil(num_images / args.num_jobs);

%2. Workout start and end indices for job
start_i	= (args.task_id-1)*job_size + 1;
end_i	= min(args.task_id*job_size, num_images);
display(['Extracting vessel centres from images ' num2str(start_i) ' to ' num2str(end_i)]);

%Loop though each image
for i_im = start_i:end_i
    im_name = pred_list(i_im).name(1:end-9);
    
    if exist([centre_dir im_name '_vc.mat'], 'file');
        display(['Skipping image ' num2str(i_im) ', ' centre_dir im_name '_vc.mat already exists. Rerun with overwrite=1 if necessary']);
        continue;
    end
    
    display(['Processing image ' num2str(i_im) ', ' datestr(now)]);  
    
    %Load in images and vessel markup
    try
        vessel_prob = u_load([prob_dir pred_list(i_im).name]);
        vessel_ori = u_load([ori_dir pred_list(i_im).name]);
        vessel_width = u_load([width_dir pred_list(i_im).name]);
        
        if ~isempty(fov_mask_dir)
            fov_mask = u_load([fov_mask_dir pred_list(i_im).name(1:end-9) '_f_mask.mat']);
            vessel_prob(~fov_mask) = 0;
            clear fov_mask;
        end
        
    catch last_err
        display(last_err.message);
        continue;
    end
    
    [vessel_centre] = extract_vessel_centres(vessel_prob, vessel_ori, vessel_width,...
        'prob_sigma',           args.prob_sigma,...
        'ori_sigma',            args.ori_sigma,...
        'width_sigma',          args.width_sigma,...
        'strong_vessel_thresh', args.strong_vessel_thresh,...
        'weak_vessel_thresh',   args.weak_vessel_thresh); %#ok
    [nrows ncols] = size(vessel_prob); %#ok
    save([centre_dir im_name '_vc.mat'], 'vessel_centre', 'nrows', 'ncols');    
   
end    
    
    
