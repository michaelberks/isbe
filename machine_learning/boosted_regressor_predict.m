function [y_fit] = boosted_regressor_predict(boostreg, X)
%BOOSTED_REGRESSOR_PREDICT Predict output using boosted regression
%   [y_fit] = boosted_regressor_predict(boostreg, X)
%
% Inputs:
%      boostreg - structure containing boosted regressor parameters:
%                 as created by BOOSTED_REGRESSOR_TRAIN
%
%      X - Input data
%
%
% Outputs:
%      y_fit - Predicted values of y
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 10-Mar-2011
% Author: Phil Tresadern 
% Email : philip.tresadern@manchester.ac.uk 
% Phone : +44 (0)161 275 5114 
% Copyright: (C) University of Manchester 

[Nx,Dx] = size(X);

% find out which samples fall in which bins
bin_edges = boostreg.bin_edges;
n_edges = size(bin_edges,1);
X_bins = zeros(size(X));
for d = 1:size(X,2)
	xx = X(:,d*ones(1,n_edges));
	be = bin_edges(:,d*ones(1,Nx))';
	X_bins(:,d) = sum(xx>be,2)+1;
end

% initialize
switch boostreg.output_type
	case 'raw',		y_fit = zeros(Nx,1);
	case 'angle',	y_fit = ones(Nx,1); % theta = 0
end

% make prediction
for L = 1:boostreg.n_levels
	input_var = boostreg.levels(L).input;
	y = boostreg.levels(L).y(X_bins(:,input_var));
	switch boostreg.output_type
		case 'raw',		y_fit = y_fit + y(:);
		case 'angle',	y_fit = y_fit .* y(:); % add angles
	end
end
