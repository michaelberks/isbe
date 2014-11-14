function [linreg] = pt_lin_reg_train(X, y, varargin)
%PT_LIN_REG_TRAIN Train a linear regressor from training pairs (X,y)
%   [X, y] = linear_regressor_train(beta, y_off)
%
% Inputs:
%      X:      N x D matrix of input data, where each row is a
%              datapoint consisting of d input variables
%
%      y:      N x d vector of class labels for each data point
%
%
%      regularization:
%              Method for regularizing regression coefficients
%              Particularly useful if training examples are limited
%
%      reg_parameter:
%              Parameter controlling behaviour of regularizer (and is
%              therefore dependent on regularizer chosen)
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
%   TREE_REG_TRAIN
%
% Created: 10-Mar-2011
% Author: Phil Tresadern 
% Email : philip.tresadern@manchester.ac.uk 
% Phone : +44 (0)161 275 5114 
% Copyright: (C) University of Manchester 

% Unpack the arguments:
args = u_packargs(varargin,... % the user's input
			 '0', ... % non-strict mode
             'regularization', 'none', ...
			 'reg_parameter', [] );

f_debug = (nargin==0 && nargout==0);
if f_debug
	X = rand(250000,108);
	y = rand(250000,1);
	y(1000,:) = nan;
end
		 
[Nx,Dx] = size(X);
[Ny,Dy] = size(y);

if (Nx~=Ny)
	error('Different number of inputs to outputs');
end
N = Nx; % = Ny

% use largescale method if big matrix
f_largescale = (numel(X)>1e7);

% find and move any NaNs in outputs
% if matrix is very big, deleting rows causes menory upsets
nan_inds = find(any(isnan(y),2));
if ~isempty(nan_inds)
	y(nan_inds,:) = y(end-length(nan_inds)+1:end,:);
	X(nan_inds,:) = X(end-length(nan_inds)+1:end,:);
	N = N-length(nan_inds);
end

% normalize outputs
y_mean = mean(y(1:N,:));
y_sd = std(y(1:N,:));
y = (y(1:N,:)-y_mean(ones(N,1),:)) ./ y_sd(ones(N,1),:);

% normalize inputs
if f_largescale
	% we do it this horrible way because std() and var() both replicate the
	% mean so that it's the same size as X - not good if X is very large
	X_mean = zeros(1,Dx);
	X_sd = zeros(1,Dx);
	for c = 1:size(X,2)
		X_mean(c) = mean(X(1:N,c));
		X_sd(c) = sum((X(1:N,c)-X_mean(c)).^2) / (N-1);
		X(1:N,c) = (X(1:N,c)-X_mean(c))/X_sd(c);
	end
else
	X_mean = mean(X(1:N,:));
	X_sd = std(X(1:N,:));
	X(1:N,:) = (X(1:N,:)-X_mean(ones(N,1),:)) ./ X_sd(ones(N,1),:);
end

% augment X with ones to account for any constant offset
switch args.regularization
	case 'none',
		if f_largescale
			rows = 1:10000;
			XtX = zeros(Dx,Dx);
			Xty = zeros(Dx,1);
			while ~isempty(rows)
				XtX = XtX + X(rows,:)'*X(rows,:);
				Xty = Xty + X(rows,:)'*y(rows);
				rows = rows+10000;
				rows(rows>N) = [];
			end
			beta = XtX \ Xty;
		else
			beta = X(1:N,:) \ y;
		end
		
	case 'pcr',
		% principal component regression
		
		% limit minimum variance proportion to range [0,1]
		min_var = min(max(args.reg_parameter,0),1);
		
		% compute SVD
		[U,S,V] = svd(X(1:N,:),'econ');
		
		% compute number of dimensions to keep
		cum_var = cumsum(diag(S));
		cum_var_prop = cum_var/cum_var(end);
		n_dims = sum(cum_var_prop<=min_var)+1;
		
		% throw away unnecessary dimensions
		if (n_dims<Dx)
			U = U(:,1:n_dims);
			S = S(1:n_dims,1:n_dims);
			V = V(:,1:n_dims);
		end
		
		% compute linear regression in n_dims dimensions
		X_reg = U*S;
		beta = X_reg \ y;
		
		% go back to Dx dimensions
		beta = V * beta;
		
	case 'ridge',
		% ridge regression
		
		% get lambda
		lambda = args.reg_parameter;
		
		% regularize (accounting for number of samples for stability)
		XtX_reg = X'*X + Nx*lambda*eye(Dx);
		
		% do linear regression
		beta = XtX_reg \ (X'*y);
		
% 	case 'lasso',
%		% lasso regression
% 
% 		% reg_parameter is proportion of t0 to use
% 		beta0 = X \ y;
% 		t0 = sum(abs(beta0))
% 		t = min(max(args.reg_parameter,0),1) * t0;
% 
% 		% precompute these products
% 		XtX = X'*X; ytX = y'*X;
% 
% 		% set up quadratic program for beta = [beta_+; beta_-]
% 		% see Tibshirani, JRSSB96, p279
% 		H	= [XtX -XtX; -XtX XtX];
% 		f	= [ytX -ytX];
% 		A	= ones(1,2*Dx);
% 		b	= t;
% 		LB	= zeros(1,2*Dx);
% 		UB	= inf(1,2*Dx);
% 		beta_aug0	= zeros(2*Dx,1);
% 		
% 		% solve quadratic program
% 		beta_aug = quadprog(H,f,A,b,[],[],LB,UB,beta_aug0);
% 		
% 		% reassemble into a single beta vector
% 		beta = reshape(beta_aug,[Dx,2]);
% 		beta = beta(:,1)-beta(:,2);
		
	otherwise,
		error(['Unknown regularization method: ',args.regularization]);
end

% compensate for scaling in X and y
beta = beta ./ X_sd(ones(Dy,1),:)';
beta = beta .* y_sd(ones(Dx,1),:);

% add constant offset to account for shifts in X and y
beta = [y_mean-X_mean*beta; beta];

linreg = struct('beta',beta,...
				'regularization',args.regularization,...
				'reg_parameter',args.reg_parameter);
			
if f_debug
	clear;
end
