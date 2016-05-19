function [] = csf_synthetic_flow(varargin)
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
    'vessel_num',		unixenv('SGE_TASK_ID',1), ...
    'data_root',        unixenv('DATA_ROOT', 'C:\isbe\nailfold\data\rsa_study\set12g'), ...
    'model_root',       unixenv('MODEL_ROOT', 'C:\isbe\nailfold\models'), ...
    'flow_map_dir',     unixenv('FLOW_MAP_DIR', 'flow_maps'), ...
    'flow_data_dir',	unixenv('FLOW_DATA_DIR', 'syn_flow_data'),...
    'flow_results_dir',	unixenv('FLOW_RESULTS_DIR', 'syn_flow_results'),...
    'flow_metrics_dir',	unixenv('FLOW_METRICS_DIR', 'syn_flow_metrics'),...
    'flow_params_dir',	unixenv('FLOW_PARAMS_DIR', 'syn_flow_params'),...
    'plot',             0);

display(args);

%Construct paths to images and predictor
flow_map_dir = [args.data_root '/' args.flow_map_dir '/'];
flow_data_dir = [args.data_root '/' args.flow_data_dir '/'];
flow_results_dir = [args.data_root '/' args.flow_results_dir '/'];
flow_metrics_dir = [args.data_root '/' args.flow_metrics_dir '/'];
flow_params_dir = [args.data_root '/' args.flow_params_dir '/'];

create_folder(flow_data_dir);
create_folder(flow_results_dir);
create_folder(flow_metrics_dir);
create_folder(flow_params_dir);

flow_list = dir([flow_map_dir '*.mat']);

vessel_name = flow_list(args.vessel_num).name(1:end-7);
display(['Processing vessel ' num2str(args.vessel_num) ': ' vessel_name]);

noise = [0.01*(2:6) 0.01*(2:6)];
flow = [3*ones(1,5) 6*ones(1,5)];

for i_rpt = 11:20

    %Randomly sample max flow and noise properties
    vessel_contrast = 1;%2*rand;
    speckle_noise = noise(i_rpt-10);%0.02*rand + 0.05;
    mean_flow = flow(i_rpt-10);%1+9*rand;

    %Make synthetic frames and save
    [frames, cell_positions] = make_synthetic_flow_frames(...
        vessel_name, mean_flow,...
        'vessel_contrast', vessel_contrast,...
        'speckle_noise', speckle_noise,...
        'flow_map_dir', flow_map_dir);        
    save([flow_data_dir vessel_name '_rpt' zerostr(i_rpt, 2)], ...
        'frames', 'cell_positions');
    save([flow_params_dir vessel_name '_rpt' zerostr(i_rpt, 2)], ...
        'vessel_contrast', 'mean_flow', 'speckle_noise');

    %Compute flow estimates  
    flow_results = [];
    [flow_results.flowPyramidEst, flow_results.flowConfidence] = ...
        estimate_flow_multilevel(255*frames, [], [], 1:4);

    save([flow_results_dir vessel_name '_rpt' zerostr(i_rpt, 2)], 'flow_results');

    %Compute flow metrics
    [flow_metrics] =  compute_vessel_flow_rf(...
        'frames', frames,...
        'flow_results', flow_results.flowPyramidEst{1},...
        'rescale_factor', 4/3,...
        'patch_contrast_scale', 40,...
        'model_root', [args.model_root '/']);
    save([flow_metrics_dir vessel_name '_rpt' zerostr(i_rpt, 2)], 'flow_metrics');
    
    if args.plot
        figure; 
        subplot(1,3,1); imgray(mean(frames,3));
        title(['MFR: true = ' num2str(mean_flow,3) ' , est = ' num2str(flow_metrics.weighted_flow_rate, 3)]);
        subplot(1,3,2); imgray(flow_metrics.vessel_pred);
        plot(flow_metrics.apex_xy(:,1), flow_metrics.apex_xy(:,2),'rx');
        title(['MWE = ' num2str(flow_metrics.mean_weighted_error, 3)]);
        subplot(1,3,3); imgray(complex2rgb(flow_results.flowPyramidEst{1}, [], [], [], 1));
        title(['Flow ratio = ' num2str(flow_metrics.vessel_flow / flow_metrics.background_flow, 3)]); 
    end
        
end
