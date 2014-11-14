function [arg_value] = load_value_from(arg_text, arg_name, default_value)
%LOAD_VALUE_FROM *Insert a one line summary here*
%   [arg_value] = load_value_from(arg_text, arg_name, default_value)
%
% Inputs:
%      arg_text - *Insert description of input variable here*
%
%      arg_name - *Insert description of input variable here*
%
%      default_value - *Insert description of input variable here*
%
%
% Outputs:
%      arg_value - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 10-Jun-2014
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

if ~exist('default_value', 'var')
    default_value = [];
end

try
    arg_value = arg_text{2}{strncmpi(arg_text{1}, arg_name, length(arg_name))};
    if isnumeric(default_value)
       arg_value = str2num(arg_value); %#ok
    end
catch %#ok
    arg_value = default_value;
end
