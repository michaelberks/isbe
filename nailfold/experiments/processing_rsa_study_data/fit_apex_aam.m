function [] = fit_apex_aam(varargin)
%ANALYSE_APEX_IMAGE_PROPERTIES *Insert a one line summary here*
%   [] = analyse_apex_image_properties()
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
    'do_initialisation',    1,...
    'do_aam_fit',           1,...
    'data_dir',             [nailfoldroot 'data/rsa_study/test/'],...
    'candidates_dir',       'apex_maps/frog/post_merged/',...
    'image_dir',            'images/',...
    'ori_dir',              'rf_regression/259076/',...
    'width_dir',            'rf_regression/257847/',...
    'aam_dir',              'aam/',...
    'aam_name',             'aam/orig/2/vessel_apex_orig.smd',...
    'model_dir',            [nailfoldroot 'data/rsa_study/models/apex_templates/'],...
    'ori_sigma',            0,...
    'width_sigma',          2,...
    'strong_vessel_thresh', 0.25,...
    'weak_vessel_thresh',   0,...
    'delete_candidate_patches', 0);

aam_dir = [args.data_dir args.aam_dir];
candidates_dir = [args.data_dir args.candidates_dir];
image_dir = [args.data_dir args.image_dir];
ori_dir = [args.data_dir 'predictions/orientation/' args.ori_dir];
width_dir = [args.data_dir 'predictions/width/' args.width_dir];

test_list = dir([candidates_dir '*.mat']);
num_images = length(test_list);

%1. Workout number of images in job
job_size = ceil(num_images / args.num_jobs);

%2. Workout start and end indices for job
start_i	= (args.task_id-1)*job_size + 1;
end_i	= min(args.task_id*job_size, num_images);
display(['Extracting vessel centres from images ' num2str(start_i) ' to ' num2str(end_i)]);

if args.do_initialisation
    %Load in model templates
    load([args.model_dir 'mean_shape.mat'], 'mean_shape');   
end

%Loop though each image
for i_im = start_i:end_i
    
    im_num = test_list(i_im).name(1:6);
    display(['Processing image ' num2str(i_im) ' of ' num2str(num_images)]); 
    
    if args.do_initialisation
        
        %Load in images and vessel markup
        load([candidates_dir im_num '_candidates.mat'], 'candidate_xy');
        vessel_ori = u_load([ori_dir im_num '_pred.mat']);
        vessel_width = u_load([width_dir im_num '_pred.mat']);
    
        image_path = [image_dir im_num '.mat'];

        initialise_aam_candidates(image_path, vessel_width, vessel_ori, candidate_xy, ...
            'aam_dir', [aam_dir im_num '/'],...
            'mean_shape', mean_shape,...
            'base_width', 15,...
            'width_sigma', 2,...
            'ori_sigma', 0,...
            'max_num_candidates', inf,...
            'thetas', linspace(-pi/8, pi/8, 20),...
            'scales', [0.8 0.9 1.0 1.1 1.2],...
            'debug', 0);
    end
    if args.do_aam_fit
    
        fit_aam_to_candidates(...
            'aam_dir', [aam_dir im_num '/'],...
            'aam_exe', 'ncm_sandpit_mb',...
            'aam_path', [args.model_dir args.aam_name],...
            'delete_candidate_patches', args.delete_candidate_patches);
    end
       
end
        