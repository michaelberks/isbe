function [] = make_detection_results_struc(varargin)
%MAKE_DETECTION_RESULTS_STRUC *Insert a one line summary here*
%   [] = make_detection_results_struc(varargin)
%
% MAKE_DETECTION_RESULTS_STRUC uses the U_PACKARGS interface function
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
% Created: 31-Oct-2013
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

% Unpack the arguments:
args = u_packargs(varargin, 0,...
    {...
    'results_name',...
    'candidates_dir'    },... %mandatory arguments
    'apex_gt_dir', [nailfoldroot 'data/rsa_study/test/apex_gt/'],...
    'results_dir', [nailfoldroot 'data/rsa_study/test/results/'],...
    'prob_dir', [nailfoldroot 'data/rsa_study/test/predictions/detection/rf_classification/257273/'],...
    'selected_gt', [],...
    'selected_candidates', [],...
    'selected_prob', [], ...
    'prob_smoothing_sigma', 2, ...
    'do_plot', 0);

clear varargin;

%Get list of ground truth and canddiates
gt_list = dir([args.apex_gt_dir '*.mat']);
if isempty(args.selected_gt)
    args.selected_gt = 1:length(gt_list);
end

candidates_list = dir([args.candidates_dir '*.mat']);
if isempty(args.selected_candidates)
    args.selected_candidates = 1:length(candidates_list);
end

if length(args.selected_candidates) ~= length(args.selected_gt)
    error('Number of selected candidates and ground truth files do not match');
end

gt_list = gt_list(args.selected_gt);
candidates_list = candidates_list(args.selected_candidates);
num_images = length(gt_list);

if ~isempty(args.prob_dir)
    prob_list = dir([args.prob_dir '*.mat']);
    if isempty(args.selected_prob)
        args.selected_prob = 1:length(prob_list);
    end
    if num_images ~= length(args.selected_prob)
        error('Number of selected candidates and vessel probability do not match');
    end
    prob_list = prob_list(args.selected_prob);
end

%
detections = cell(num_images, 3);
%
for i_im = 1:num_images;
    
    display(['Processing image ' num2str(i_im) ' of ' num2str(num_images)]);
    
    
    load([args.apex_gt_dir gt_list(i_im).name]);
    load([args.candidates_dir candidates_list(i_im).name],...
        'candidate_xy', 'candidate_scores');
       
    if ~isempty(args.prob_dir)
        vessel_prob = u_load([args.prob_dir prob_list(i_im).name]);
        
        if args.prob_smoothing_sigma
            g = gaussian_filters_1d(args.prob_smoothing_sigma);
            g = g / sum(g);
            vessel_prob = conv2(g', g, vessel_prob, 'same');
        end
        
    else
        vessel_prob = [];
    end
    
    do_plot = args.do_plot && i_im <= 20;
    
    [detections{i_im,1} detections{i_im,2} detections{i_im,3}] =...
        evaluate_apex_candidates(apex_xy, candidate_xy, apex_widths, vessel_prob, ...
        [], [], do_plot);    
    
end
%
if isempty(args.results_name)
    args.results_name = ['detection_results_' datestr(now,30)];
end

create_folder(args.results_dir);
save([args.results_dir args.results_name '.mat'], ...
    'detections', 'gt_list', 'candidates_list');