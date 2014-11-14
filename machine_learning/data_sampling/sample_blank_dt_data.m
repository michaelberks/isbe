function [training_data training_labels] = sample_blank_dt_data(sampling_args)
%SAMPLE_DT_DATA_TEMP *Insert a one line summary here*
%   [] = sample_dt_data_temp(varargin)
%
% SAMPLE_DT_DATA_TEMP uses the U_PACKARGS interface function
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
% Created: 03-Nov-2010
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 


num_blank_samples = round(0.05*sampling_args.num_samples);

sampling_args.num_samples = sampling_args.num_samples - num_blank_samples;

[training_data training_labels] = generate_line_training_data(sampling_args);
    
if ispc
    blank_bg_dir = 'blank_bgs';
else
    [z blank_bg_dir] = unix('echo $BLANK_BG_DIR'); blank_bg_dir(end) = [];
end
blank_bg_dir = [asymmetryroot 'data/synthetic_backgrounds/' blank_bg_dir '/'];

num_blank_bgs = round(num_blank_samples / sampling_args.pts_per_image);
blank_bg_list = dir([blank_bg_dir '*.mat']);
num_blank_bgs = min(num_blank_bgs, length(blank_bg_list));
rp = randperm(length(blank_bg_list));
blank_bg_list = blank_bg_list(rp(1:num_blank_bgs));

sampling_args2.num_samples = num_blank_samples;
sampling_args2.image_dir = blank_bg_dir;
sampling_args2.image_list = blank_bg_list;
sampling_args2.win_size = sampling_args.win_size;
sampling_args2.num_levels = sampling_args.num_levels;
sampling_args2.feature_type = sampling_args.feature_type;
sampling_args2.feature_shape = sampling_args.feature_shape;
sampling_args2.rotate = sampling_args.rotate;
sampling_args2.do_max = sampling_args.do_max;
sampling_args2.use_nag = sampling_args.use_nag;

training_data = [training_data; ...
                 sample_image_dt_data(sampling_args2)];
training_labels = [training_labels;...
                   false(num_blank_samples, 1)];

