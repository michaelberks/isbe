function [apex_measures, apex_data] = make_apex_measures_from_cxx(data_dir, image_name, varargin)
%MAKE_APEX_MEASURES_FROM_CXX *Insert a one line summary here*
%   [apex_measures] = make_apex_measures_from_cxx(cxx_apex_data_file)
%
% Inputs:
%      cxx_apex_data_file - *Insert description of input variable here*
%
%
% Outputs:
%      apex_measures - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 18-May-2016
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester

args = u_packargs(varargin,... % the user's input
    '0', ... % non-strict mode
    'capillary_type',               {'distal', 'nondistal'},...
    'vessel_prediction_dir',        'predictions/vessels',...
    'orientation_prediction_dir',   'predictions/orientations',...
    'width_prediction_dir',         'predictions/widths',...
    'apex_prediction_dir',          'predictions/apexes',...
    'vessel_prediction_ext',        '_pred.png',...
    'orientation_prediction_ext',   '_pred.png',...
    'width_prediction_ext',         '_pred.png',...
    'apex_candidates_ext',          '_apex_candidates.txt',...
    'prob_sigma',                   2,...
	'ori_sigma',                    0,...
	'width_sigma',                  2,...
	'num_c_pts',                    20,...
	'max_dist',                     120,...
	'border_sz',                    16,...
	'num_ori_bins',                 36,...
	'connect_thresh',               0.5,...
	'base_width',                   20,...
	'width_predictor',              [],...
	'plot',                         0);


%Set up output containers
apex_measures.distal = [];
apex_measures.nondistal = [];
apex_data = [];

%Set up args for extracting apex data
apex_args.prob_sigma = args.prob_sigma;
apex_args.ori_sigma = args.ori_sigma;
apex_args.width_sigma = args.width_sigma;
apex_args.num_c_pts = args.num_c_pts;
apex_args.max_dist = args.max_dist;
apex_args.border_sz = args.border_sz;
apex_args.num_ori_bins = args.num_ori_bins;
apex_args.connect_thresh = args.connect_thresh;
apex_args.base_width = args.base_width;
apex_args.width_predictor = args.width_predictor;
apex_args.plot = args.plot > 1;

candidates_filename = [data_dir '/' args.apex_prediction_dir '/'...
    image_name args.apex_candidates_ext];
if ~exist(candidates_filename, 'file')
    display(['Missing apex file ' candidates_filename]);
    return;
end
    
apex_data = load(candidates_filename); 
if isempty(apex_data)
    display([candidates_filename ' contains no apex data']);
    return;
end
        
vessel_prob_name = [data_dir '/' args.vessel_prediction_dir '/'...
    image_name args.vessel_prediction_ext];
vessel_ori_name = [data_dir '/' args.orientation_prediction_dir '/'...
    image_name args.orientation_prediction_ext];
vessel_width_name = [data_dir '/' args.width_prediction_dir '/'...
    image_name args.width_prediction_ext];

if ~exist(vessel_prob_name, 'file')
    display([vessel_prob_name ' not found']);
    return;
elseif ~exist(vessel_ori_name, 'file')
    display([vessel_ori_name ' not found']);
    return;
elseif ~exist(vessel_width_name, 'file')
    display([vessel_width_name ' not found']);
    return;
end

display(['processing sequence ' candidates_filename]);

vessel_prob = imread(vessel_prob_name);
vessel_prob = double(vessel_prob)/100;
vessel_ori = imread(vessel_ori_name);
vessel_ori = rgb2complex(vessel_ori, [], 1, [], 0);
vessel_width = imread(vessel_width_name);
vessel_width = double(vessel_width);   

candidate_xy = apex_data(:,1:2);
selected_distal = apex_data(:,11) > 0;
selected_nondistal = apex_data(:,12) > 0;

for i_type = 1:length(args.capillary_type)

    if strcmpi(args.capillary_type{i_type}, 'distal')
        selected_idx = selected_distal;
    elseif strcmpi(args.capillary_type{i_type}, 'nondistal')
        selected_idx = selected_nondistal;       
    else
        display(['Capillary type ' args.capillary_type{i_type} ' not recognised']);
        continue;
    end

    [measures_struc] = extract_apex_measures_newest(...
        vessel_prob, vessel_ori, vessel_width, vessel_prob, candidate_xy(selected_idx,:),... 
        apex_args);

    %Copying the fields in this way means we don't overwrite existing
    %data
    fnames = fieldnames(measures_struc);
    for i_f = 1:length(fnames)
        apex_measures.(args.capillary_type{i_type}).(fnames{i_f}) = measures_struc.(fnames{i_f});
    end

    apex_measures.(args.capillary_type{i_type}).candidate_scores = apex_data(selected_idx,7);
    apex_measures.(args.capillary_type{i_type}).candidate_displacements = apex_data(selected_idx,8);
    apex_measures.(args.capillary_type{i_type}).candidate_class_probs = apex_data(selected_idx,10);

end

if args.plot
    figure;
    subplot(2,1,1); imgray(vessel_prob);
    plot(candidate_xy(selected_distal,1),...
         candidate_xy(selected_distal,2), 'rx', 'markersize', 10);
    plot(candidate_xy(selected_nondistal,1),...
         candidate_xy(selected_nondistal,2), 'gx');
    plot(candidate_xy(~selected_distal & ~selected_nondistal,1),...
         candidate_xy(~selected_distal & ~selected_nondistal,2), 'c.', 'markersize', 4);
    
    subplot(2,1,2); imgray(complex2rgb(vessel_ori));
end
