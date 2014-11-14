function [y_fit] = linear_regressor_predict(linreg, X)
%LINEAR_REGRESSOR_PREDICT Predict output using linear regression
%   [y_fit] = linear_regressor_predict(beta, y_off, X)
%
% Inputs:
%      linreg - structure containing linear regressor parameters:
%                 beta - D x d matrix of regression coefficients
%                 y_off - 1 x d vector of output offsets
%      as created by LINEAR_REGRESSOR_TRAIN
%
%      X - Input data
%
%
% Outputs:
%      y_fit - Predicted values of y (= y_off + X*beta)
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
if (size(linreg.beta,1) ~= (Dx+1))
	error('Wrong number of input parameters');
end

y_fit = [ones(Nx,1) X]*linreg.beta;
	
