function [exp_sample] = exp_rand_sample(mu, sample_size)
%EXP_RAND_SAMPLE *Insert a one line summary here*
%   [exp_sample] = exp_rand_sample(mu,sample_size)
%
% Inputs:
%      mu- *Insert description of input variable here*
%
%      sample_size- *Insert description of input variable here*
%
%
% Outputs:
%      exp_sample- *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 16-Apr-2010
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester
if nargin < 2
    sample_size = 1;
end
uni_sample = rand(sample_size);
exp_sample = -mu*log(1-uni_sample);
