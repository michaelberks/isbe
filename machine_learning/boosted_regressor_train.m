function [predictor] = boosted_regressor_train(varargin)
%boosted_REGRESSOR_TRAIN Sample data and build a boosted regressor
%   [X, y] = boosted_regressor_train(beta, y_off)
%
% Inputs:
%   sampling_method: Function name used to extract each set of training
%       data from the global population. It is assumed this produces an N x D 
%       matrix of training data (where each row is a datapoint consisting of D
%       input variables) and an N x 1 vector of class labels for each data
%       point
%
%   sampling_method_args: structure of arguments used by the sampling
%       method
%
%   tree_dir: path to the directory to which regressor should be saved
%
%   tree_root: root folder for regressor
%
% Outputs:
%      boostreg - structure containing boosted regressor parameters:
%                 beta - D x d matrix of regression coefficients
%                 y_off - 1 x d vector of output offsets
%
%
% Example:
%
% Notes:
%
% See also:
%   RANDOM_FOREST_REG_TRAIN
%
% Created: 10-Mar-2011
% Author: Phil Tresadern 
% Email : philip.tresadern@manchester.ac.uk 
% Phone : +44 (0)161 275 5114 
% Copyright: (C) University of Manchester 

% Unpack the arguments:
args = u_packargs(varargin,... % the user's input
    '0', ... % non-strict mode
    ... % Mandatory arguments
   {'sampling_args',... % for sampling training data
    'decomposition_args',... % for sampling training features
    'predictor_args'}, ...
    ... % Optional arguments
    'tree_root', []);

%Get input arguments for whole forest from fields of args
sampling_args = args.sampling_args;
sampling_method = sampling_args.sampling_method;

%Compute individual classification trees for boostrap samples of X
display('Building boosted regressor');

% Get training data
sampling_method_args.task_id = 1;
tic;
[X y] = feval(sampling_method, args);
t = toc;
display(['Time sampling data = ', num2str(t)]);

% build a boosted regressor for the sampled data set
tic;
predictor = pt_boosted_reg_train(X, y, args.predictor_args); %#ok
t = toc;
display(['Time building boosted regressor = ', num2str(t)]);
clear X y;

% save fields to predictor output structure
predictor_args = args.predictor_args;
predictor.args = predictor_args;

% Set filename at which to save tree
create_folder(predictor_args.model_dir);
save(predictor_args.predictor_path, 'predictor');
