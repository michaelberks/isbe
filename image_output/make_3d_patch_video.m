function [] = make_3d_patch_video(input_patch, gif_name, varargin)
%MAKE_3D_PATCH_VIDEO *Insert a one line summary here*
%   [] = make_3d_patch_video(varargin)
%
% MAKE_3D_PATCH_VIDEO uses the U_PACKARGS interface function
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
% Created: 14-May-2013
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

% Unpack the arguments:
args = u_packargs(varargin, '0', ...
    'num_frames', 72,...
    'fig_bg_color', 'white',...
    'colormap', gray(256),...
    'edgecolor', 'black',...
    'elevation', 45,...
    'camera_view_angle', 8,...
    'delay_time', 0.05);
clear varargin;

angle_step = 360 / args.num_frames;
azimuth_angles = 0:angle_step:(args.num_frames-1)*angle_step;


f1 = figure(...
    'windowstyle', 'normal', 'color', args.fig_bg_color,...
    'position', [50 50 800 800]); 
surf(input_patch, 'edgecolor', 'k'); 
colormap(args.colormap);

for az = azimuth_angles
    
    view(az, args.elevation);
    set(gca, 'CameraViewAngle', args.camera_view_angle);
    axis off;
    
    frame1 = getframe(f1);
    gif1 = frame2im(frame1);
    [gif1a map] = rgb2ind(gif1, 2^16);
    
    if ~az
        imwrite(gif1a, map, gif_name, 'gif', 'WriteMode', 'overwrite', 'Loop', Inf, 'DelayTime', args.delay_time);
    else
        imwrite(gif1a, map, gif_name, 'gif', 'WriteMode', 'append', 'DelayTime', args.delay_time);
    end
end