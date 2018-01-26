function [data_out] = add_rician_noise(data_in, sigma)
%ADD_RICIAN_NOISE *Insert a one line summary here*
%   [data_out] = add_rician_noise(data_in, sigma)
%
% Inputs:
%      data_in - *Insert description of input variable here*
%
%      sigma - *Insert description of input variable here*
%
%
% Outputs:
%      data_out - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 25-Jul-2017
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

data_out = abs(data_in +...
    sigma.*randn(size(data_in)) +...
    sqrt(-1).*sigma.*randn(size(data_in)) );
