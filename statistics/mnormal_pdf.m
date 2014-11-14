function [p] = mnormal_pdf(X, mu, sigma)
%MNORMAL_PDF *Insert a one line summary here*
%   [p] = mnormal_pdf(X, mu, sigma)
%
% Inputs:
%      X - *Insert description of input variable here*
%
%      mu - *Insert description of input variable here*
%
%      sigma - *Insert description of input variable here*
%
%
% Outputs:
%      p - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 04-Apr-2013
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester
k = size(sigma,1);

A = 1 / ( (2*pi)^(k/2) * sqrt(det(sigma)) );
p = A*exp(-0.5*(X-mu)/sigma*(X-mu)' );
