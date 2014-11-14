function [predictor] = linear_regressor_train(varargin)
%LINEAR_REGRESSOR_TRAIN Sample data and build a linear regressor
%   [X, y] = linear_regressor_train(beta, y_off)
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
%      linreg - structure containing linear regressor parameters:
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
    'nonstrict', ...
    ... % mandatory arguments
   {'sampling_args', ...
    'predictor_args', ...
    'decomposition_args'});

%Get input arguments for whole forest from fields of args
sampling_args = args.sampling_args;
sampling_method = sampling_args.sampling_method;

%Compute individual classification trees for boostrap samples of X
display('Building linear regressor');

% Get training data
sampling_args.task_id = 1;
tic;
[X y] = feval(sampling_method, args);
t = toc;
display(['Time sampling data = ', num2str(t)]);

if ispc && strcmp(get_username,'ptresadern')
% 	save([linreg_dir 'traindata0001.mat'], 'X','y');
end

% build a linear regressor for the sampled data set
tic;
if any(~isreal(y))
    % standard deviations aren't well defined for complex numbers so do
    % real and imaginary parts independently. (This probably is not correct
    % for all applications, but will do for our purposes.)
    predictor = pt_lin_reg_train(X, real(y));
    predictor_imag = pt_lin_reg_train(X, imag(y));
    predictor.beta = complex(predictor.beta, predictor_imag.beta);
else
    predictor = pt_lin_reg_train(X, y); %#ok
end
t = toc;
display(['Time building linear regressor = ', num2str(t)]);
clear X y;

% save fields to predictor output structure
predictor_args = args.predictor_args;
predictor.args = predictor_args;

% Set filename at which to save tree
create_folder(predictor_args.model_dir);
save(predictor_args.predictor_path, 'predictor');
