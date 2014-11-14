function [training_data training_labels] = sample_saved_training_data(varargin)
%sample_training_data *Insert a one line summary here*
%   [] = sample_training_data(varargin)
%
% sample_training_data uses the U_PACKARGS interface function
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
% Created: 17-Aug-2010
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester

f_debug = (nargin == 0 && nargout == 0);
if f_debug, test_script(); return; end

[training_data training_labels] = func(varargin{:});


%% The function
function [training_data training_labels] = func(varargin)
args = u_packargs(varargin,... % the user's input
    '0', ... % non-strict mode
    ... % Mandatory arguments
   {'sampling_args'},...
    'quiet', true);
clear varargin;

% Workout feature length from decomposition parameters
[training_data training_labels boot_idx] = bootstrap_train_data(...
    args.sampling_args.X, args.sampling_args.y);

% All training features/labels are now in. Save the parameter lists if they
% exist - a bit messy to have this here, but to keep the tree training code
% general the generate training data function can only return X & y. So we
% use the model wrapper function to pass a path that the model will be
% saved to, so the parameters/sampling points data can also be saved to
% this path
if isfield(args.sampling_args, 'sampled_data_dir') && ~isempty(args.sampling_args.sampled_data_dir)
    create_folder(args.sampling_args.sampled_data_dir);
    par_list = dir([args.sampling_args.sampled_data_dir '/sampled_pts*.mat']);
    tree_num = length(par_list) + 1;
    if exist('sampled_pts','var') && ~isempty(sampled_pts)
        save([args.sampling_args.sampled_data_dir '/sampled_pts' zerostr(tree_num, 3) '.mat'], 'boot_idx');
    end
end


%% Test script
function test_script()
clc;