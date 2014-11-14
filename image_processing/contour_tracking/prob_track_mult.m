function [path_map] = prob_track_mult(I_prob, I_ori, varargin)
%PROB_TRACK_MULT *Insert a one line summary here*
%   [path_map] = prob_track_mult(varargin)
%
% PROB_TRACK_MULT uses the U_PACKARGS interface function
% and so arguments to the function may be passed as name-value pairs
% or as a struct with fields with corresponding names. The field names
% are defined as:
%
% Mandatory Arguments:
%
% Optional Arguments:
%
% Outputs:
%      path_map - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 23-Feb-2012
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

% Unpack the arguments:
args = u_packargs(varargin, '0', ...
    'num_paths', 1e5,...
    'junction_map', [],...
    'step_length', 1,...
    'double_angle', 1,...
    'ignore_mask', [],...
    'plot', 0);
clear varargin;

if isempty(args.ignore_mask)
    args.ignore_mask = true(size(I_prob));
end

%Compute direction and dispersion from orientation
I_ori_theta = angle(I_ori);
I_ori_D = abs(I_ori);

%Check if we need to halve the angle
if args.double_angle
    I_ori_theta = I_ori_theta / 2;
end

%preallocate path map
path_map = zeros(size(I_prob));

for ii = 1:args.num_paths
    %Choose a starting point - possibly provide options on this?
    [y_pot x_pot] = find(I_prob > rand & args.ignore_mask);
    r_idx = ceil(rand*length(x_pot));
    x1 = x_pot(r_idx);
    y1 = y_pot(r_idx);
    
    %Compute path from that point
    [path_map_i] = prob_track(I_prob, I_ori_theta, I_ori_D, x1, y1, args.step_length, args.junction_map);
    
    %Add path to full map
    path_map = path_map + path_map_i;
end

if args.plot
    figure; imgray(log(path_map));
end
