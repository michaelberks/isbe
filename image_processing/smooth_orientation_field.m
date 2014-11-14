function [] = smooth_orientation_field(ori_map, varargin)
%SMOOTH_ORIENTATION_FIELD *Insert a one line summary here*
%   [] = smooth_orientation_field(varargin)
%
% SMOOTH_ORIENTATION_FIELD uses the U_PACKARGS interface function
% and so arguments to the function may be passed as name-value pairs
% or as a struct with fields with corresponding names. The field names
% are defined as:
%
% Mandatory Arguments:
%
% Optional Arguments:
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 01-Apr-2011
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

% Unpack the arguments:
args = u_packargs(varargin, '0', ...
    'mask', [],...
    'method', 0,...
    'sigma', 1);
clear varargin;