function [logreg] = logistic_regressor_train(varargin)
%LOGISTIC_REGRESSOR_TRAIN Sample data and build a logistic regressor
%   [X, y] = logistic_regressor_train(beta, y_off)
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
%      logreg - structure containing logistic regressor parameters:
%                 beta - D x d matrix of regression coefficients
%
%
% Example:
%
% Notes:
%
% See also:
%   RANDOM_FOREST_REG_TRAIN, PT_LIN_REG_TRAIN
%
% Created: 10-Mar-2011
% Author: Phil Tresadern 
% Email : philip.tresadern@manchester.ac.uk 
% Phone : +44 (0)161 275 5114 
% Copyright: (C) University of Manchester 

% Unpack the arguments:
args = u_packargs(varargin,... % the user's input
			 '0', ... % non-strict mode
			 {'sampling_method',... % the mandatory arguments
             'sampling_method_args',...
             'tree_dir'}, ...
             'tree_root', []);

%Get input arguments for whole forest from fields of args
sampling_method = args.sampling_method;
sampling_method_args = args.sampling_method_args;

logreg_dir = [args.tree_root args.tree_dir];
if ~isempty(logreg_dir)
    if ~strcmp(logreg_dir(end), '/') && ~strcmp(logreg_dir(end), '\');
        logreg_dir = [logreg_dir filesep];
    end
    if ~exist(logreg_dir, 'dir')
        mkdir(logreg_dir);
		if ~ispc
			fileattrib(logreg_dir,'+w','g'); % make tree dir accessible...
			fileattrib([logreg_dir,'/..'],'+w','g'); % ...and its parent, too
		end
    end
end

%Compute individual classification trees for boostrap samples of X
display('Building logistic regressor');

%% Get training data
sampling_method_args.task_id = 1;
tic;
% sample foreground pixels
sampling_method_args.data_class = 'fg';
[Xp yp] = feval(sampling_method, sampling_method_args);
% sample background pixels
sampling_method_args.data_class = 'bg';
[Xn yn] = feval(sampling_method, sampling_method_args);
% replace orientations with binary class labels
yp(:) = 1; yn(:) = 0;
% stack fg and bg data
X = [Xp; Xn]; y = [yp; yn];
t = toc;
display(['Time sampling data = ', num2str(t)]);

if ispc && strcmp(get_username,'ptresadern')
	save([logreg_dir 'traindata0001.mat'], 'X','y');
end

%% build a logistic regressor for the sampled data set
tic;
logreg = pt_log_reg_train(X, y); %#ok
t = toc;
display(['Time building logistic regressor = ', num2str(t)]);
clear X y;

% save fields to random_forest output structure
logreg.tree_dir = args.tree_dir;
logreg.tree_root = args.tree_root;
logreg.regression_method = 'logistic';

%Set filename at which to save tree
logreg_name = 'random_forest.mat';
save([logreg_dir logreg_name], 'logreg');
