function [] = csf_compute_flow(varargin)
%PREDICT_IMAGE_SET wrapper function to predict ouput for a set of images in
%a given directory
%
%
% Inputs:
%
%
% Outputs:
%
% Example: for use on hydra
%   NUM_JOBS=20 JOB_ID="'191658'" TEST_IMAGE_DIR="'lines512'" qsub -N class_im -t 1-20 -V matlab_code/trunk/hydra/predict_image_set.sh
%
% Notes:
%
% See also:
%
% Created: 03-Feb-2010
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester

% this now uses the utils/unixenv function that takes the key name of the
% environment variable and returns its corresponding value, returning the
% default value if the environment variable does not exist
%   e.g. envval = unixenv('ENVVARNAME',default_val);

warning('off', 'ASYM:unexpectedArgument');

args = u_packargs(varargin,... % the user's input
    '0', ... % non-strict mode
    'subject_id',		unixenv('SGE_TASK_ID',1), ...
    'num_jobs',			unixenv('NUM_JOBS',1), ...
    'data_root',        unixenv('DATA_ROOT', 'N:\Nailfold Capillaroscopy\wellcome'), ...
    'model_root',       unixenv('MODEL_ROOT', 'C:\isbe\nailfold\models'), ...
    'flow_data_dir',	unixenv('FLOW_DATA_DIR', 'flow_data'), ...
    'flow_results_dir',	unixenv('FLOW_RESULTS_DIR', 'flow_results'),...
    'flow_metrics_dir',	unixenv('FLOW_METRICS_DIR', 'flow_metrics'));

display(args);

%Construct paths to images and predictor
flow_data_dir = [args.data_root '/' args.flow_data_dir '/'];
flow_results_dir = [args.data_root '/' args.flow_results_dir '/'];
flow_metrics_dir = [args.data_root '/' args.flow_metrics_dir '/'];

subject_name = [zerostr(args.subject_id,3) 'wellcome'];
%flow_list = dir([flow_data_dir subject_name '*']);
flow_list = dir([flow_results_dir subject_name '*']);

num_vessels = length(flow_list);
display(['Computing flow for ' subject_name ',' num2str(num_vessels) ' vessels to process']);

for i_ves = 1:num_vessels
    display(['Processing vessel ' num2str(i_ves)]);
    
    frames = load([flow_data_dir flow_list(i_ves).name], 'cropped_frames',...
        'x_min', 'y_min', 'x_max', 'y_max', 'edge_mask', 'vessel_transforms');
    if ~isfield(frames, 'cropped_frames') || ...
            isempty(frames.cropped_frames) || ...
            ~any(frames.cropped_frames(:)) || any(size(frames.cropped_frames) < 32)
        continue;
    end
       
    %Contrast normalise the frames
    g_lims = prctile(double(frames.cropped_frames(:)), [1 99]);
    g_range = g_lims(2) - g_lims(1);
    
    %Compute flow results
    flow_results = [];
    [flow_results.flowPyramidEst, flow_results.flowConfidence] = ...
        estimate_flow_multilevel(255*(frames.cropped_frames-g_lims(1))/g_range, [], [], 1:4);   
    flow_results.g_lims = g_lims;   
    save([flow_results_dir flow_list(i_ves).name], 'flow_results');

    %Compute flow metrics
    [flow_metrics] =  compute_vessel_flow_rf(...
        'vessel_name', flow_list(i_ves).name(1:end-4),...
        'flow_data_dir',	flow_data_dir, ...
        'flow_results_dir',	flow_results_dir,...
        'model_root', [args.model_root '/']); %#ok
    save([flow_metrics_dir flow_list(i_ves).name], 'flow_metrics');


end