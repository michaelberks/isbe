function [] = quiver_orientations(ori_map, varargin)
%QUIVER_ORIENTATIONS *Insert a one line summary here*
%   [] = quiver_orientations(varargin)
%
% QUIVER_ORIENTATIONS uses the U_PACKARGS interface function
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
% Created: 17-May-2011
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

% Unpack the arguments:
args = u_packargs(varargin, '0', ...
    'line_map', [],...
    'spacing', 1,...
    'num_angles', 36,...
    'quiver_length', 4,...
    'axis', []);
clear varargin;

if isempty(args.axis)
    aa = gca;
else
    aa = args.axis;
end

if isempty(args.line_map)
    line_map = true(size(ori_map));
else
    line_map = args.line_map;
end

ori_map = mod(ori_map, pi);
ang_res = pi / args.num_angles;
arrow_colors = hsv(args.num_angles);

spacing_mask = false(size(line_map));
spacing_mask(1:args.spacing:end, 1:args.spacing:end) = true;


for ii = 0:args.num_angles
    theta = (ii - 0.5)*ang_res;
    %Get line_map of pixels that have orientation within theta range
    angle_mask = line_map & ...
                 spacing_mask & ...
                 ori_map > (theta - 0.5*ang_res) &...
                 ori_map <= (theta + 0.5*ang_res);
    [yy xx] = find(angle_mask);
    uu = args.quiver_length*args.spacing*cos(ori_map(angle_mask));
    vv = -args.quiver_length*args.spacing*sin(ori_map(angle_mask));

    quiver(xx - uu/2, yy - vv/2, uu, vv,...
        'Parent', aa,...
        'Color', arrow_colors(mod(ii, args.num_angles)+1,:),...
        'HitTest', 'off',...
        'ShowArrowHead', 'off',...
        'Autoscale', 'off');
end
