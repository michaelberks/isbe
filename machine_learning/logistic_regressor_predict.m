function [y_fit] = logistic_regressor_predict(logreg, X)
%LOGISTIC_REGRESSOR_PREDICT Predict output using logistic regression
%   [y_fit] = logistic_regressor_predict(logreg, X)
%
% Inputs:
%      logreg - structure containing logistic regressor parameters:
%                 beta - D x d matrix of regression coefficients
%				as created by LOGISTIC_REGRESSOR_TRAIN
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
if (size(logreg.beta,1) ~= (Dx+1))
	error('Wrong number of input parameters');
end

q = exp([ones(Nx,1) X]*logreg.beta);
y_fit = q./(1+q);
	
