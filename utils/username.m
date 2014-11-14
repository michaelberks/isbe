function [uname] = username()
%GET_USERNAME Find username of current user, regardless of OS
%   [username] = get_username()
%
% Inputs:
%
% Outputs:
%      username - Your username
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

if ispc
	uname = getenv('username');
elseif isunix
	uname = unixenv('USER');
end

	