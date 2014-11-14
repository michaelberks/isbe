function [env_val] = unixenv(env_key, default_val)
%UNIXENV Look up a Unix-style environment variable and return it in the
%correct format.
%   [env_val] = unixenv(env_key, default_val)
%
% Inputs:
%      env_key - Key name of the environment variable (e.g. USER)
%
%      default_val - Default value if environment variable doesn't exist
%
%
% Outputs:
%      env_val - Return value
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 03-Feb-2011
% Author: Phil Tresadern 
% Email : philip.tresadern@manchester.ac.uk 
% Phone : +44 (0)161 275 5114 
% Copyright: (C) University of Manchester 

% if no default given then assume empty
if nargin<2, default_val = []; end

% if not unix then return default value
if ispc
	env_val = default_val;
	return; 
end
	
% remove any leading $'s from key name
while ~isempty(env_key) && env_key(1)=='$'
	env_key = env_key(2:end);
end

% if env_key is empty, return empty
if isempty(env_key)
	env_val = [];
	return; 
end

% look for environment variable
[tmp env_val] = unix(['echo $',env_key]);

% remove any leading or trailing whitespace (e.g. newlines)
env_val = strtrim(env_val);

% if it doesn't exist, return default value
if isempty(env_val)
	env_val = default_val;
	return
end

if isempty(default_val)
	% determine format automatically
	
	% try converting to number
	if ~isempty(str2num(env_val))
		env_val = str2num(env_val);
    else
        %try converting to a cell
        env_val = convert2cell(env_val);
    end
else
	% use same format as default value
	
	if isnumeric(default_val) || islogical(default_val)
		% return conversion to number
		env_val = str2num(env_val);
	else
        %try converting to a cell
        env_val = convert2cell(env_val);
    end
end


