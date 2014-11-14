function [logreg] = pt_lin_reg_train(X, y, varargin)
%PT_LIN_REG_TRAIN Train a logistic regressor from training pairs (X,y)
%   [X, y] = pt_lin_reg_train(X, y, regularization, reg_parameter)
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
%      logreg - structure containing logistic regressor parameters:
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
		 
% if run as a script then do debugging
f_debug = 0;
if (nargin==0)
	f_debug = 1;
	overlap = 0.2;
	X = 2*rand(100,1);
	y = (X>1);
	X(y) = X(y)-overlap/2;
	X(~y) = X(~y) + overlap/2;
end

[Nx,Dx] = size(X);
[Ny,Dy] = size(y);

if (Nx~=Ny)
	error('Different number of inputs to outputs');
end

% add bias term to inputs
X = [ones(Nx,1) X];

% we will solve this using iteratively reweighted least squares
% see Hastie, Tibshirani & Friedman, "The Elements of Statistical
% Learning", 1st edition, pp. 95-99

% this uses an implementation for the two class case that avoids having a
% matrix W at all

% intialize regression coefficients
beta = zeros(Dx+1,1); 
% W = sparse(Nx,Nx); 

% train regressor
switch args.regularization
	case 'none',
		imax = 100;
		for i = 1:imax
			% precompute this
			Xbeta = X*beta;
			
			% estimated probabilities
			q = exp(X*beta); 
			p = q ./ (1+q); 
			
			% weight matrix
% 			W(1:Nx+1:end) = p.*(1-p);
			Wdiag = p.*(1-p);
			
			% compute the adjusted response
% 			z = Xbeta + W\(y-p);
			z = Xbeta + (y-p)./Wdiag;
			
			% compute new beta
% 			XtW = X'*W;
			XtW = zeros(size(X'));
			for j = 1:size(X,1), XtW(:,j) = Wdiag(j)*X(j,:)'; end
			beta = (XtW*X) \ (XtW*z);
			
			% compute log likelihood and stop if converged
			loglhood(i) = sum(y.*Xbeta - log(1+q));
			if (i>1),
				dl = -(loglhood(i)-loglhood(i-1))/loglhood(i-1);
				if dl<1e-3, break; end
			end
			
			if f_debug				
				figure(1); clf; 
				subplot(2,1,1); cla; hold on;
					plot(X(:,2),y,'bo');
					plot(X(:,2),p,'r.');
					ylim([0,1] + 0.1*[-1 1]);
				subplot(2,1,2); hold on;
					plot(1:i,loglhood,'b.-');
					axis([0,imax,loglhood(1),0]);
				
				pause(0.05);
			end
		end
		
	otherwise,
		error(['Unknown regularization method: ',args.regularization]);
end

logreg = struct('beta',beta,...
				'regularization',args.regularization,...
				'reg_parameter',args.reg_parameter);
