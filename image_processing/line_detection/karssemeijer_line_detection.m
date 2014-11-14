function [line_orientation, line_map, line_scale, line_strength] = karssemeijer_line_detection(image_in, varargin)
%KARSSEMEIJER_LINE_DETECTION *Insert a one line summary here*
%   [line_strength, line_orientation, line_map] = karssemeijer_line_detection(varargin)
%
% KARSSEMEIJER_LINE_DETECTION uses the U_PACKARGS interface function
% and so arguments to the function may be passed as name-value pairs
% or as a struct with fields with corresponding names. The field names
% are defined as:
%
% Mandatory Arguments:
%
% Optional Arguments:
%
% Outputs:
%      line_strength - *Insert description of input variable here*
%
%      line_orientation - *Insert description of input variable here*
%
%      line_map - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 10-Jan-2011
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

% Unpack the arguments:
args = u_packargs(varargin, '0', ...
    'line_scales', [1 2 4 8],...
    'grad_scale', 10,...
    'grad_ori_thresh', pi/6,...
    'grad_strength_thresh', 25,...
    'line_strength_thresh', 0,...
    'binary_map', 0,...
    'degrees', 1);
clear varargin;

%Compute line strength and orientation using Gaussian 2nd derivative filters
% at a range of scales
[line_strength, line_orientation, line_scale] = gaussian_2nd_derivative_line(image_in, args.line_scales);

%Compute gradient scale and orientation at a fixed coarse scale
[grad_strength, grad_orientation] = gaussian_1st_derivative_gradient(image_in, args.grad_scale);

%Discard pixels with strong gradient and similar gradient and line
%orientations as these are likely to be edges not lines. Also discard
%pixels with +ve line strength as these are negative lines (i.e. dark on
%light)
discard_map = ...
    ((abs(mb_mod(line_orientation - grad_orientation + pi/2, pi)) < args.grad_ori_thresh) &...
    (grad_strength > args.grad_strength_thresh)) | ...
    (line_strength > args.line_strength_thresh);

if args.binary_map
    %For a binary line map as output, we just invert the discard map
    line_map = ~discard_map;
else
    %Copy the absolute values of the line srength into the line map
    line_map = abs(line_strength);

    %Discard the pixels as outlined above
    line_map(discard_map) = 0;

    %Scale the line map between 1 and 0
    line_map = line_map ./ max(line_map(:));
end

line_orientation = mod(line_orientation, pi);

%if requested, covert orientations from radians to degrees
if args.degrees
    line_orientation = 180*line_orientation/pi;
end