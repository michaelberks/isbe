function [] = make_improfile_video(input_patch, gif_name, xc, yc, varargin)
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
    'profile_radius', 20,...
    'num_frames', 72,...
    'fig_bg_color', 'white',...
    'colormap', gray(256),...
    'axes_color', 'r',...
    'profile_color', 'r',...
    'delay_time', 0.05);
clear varargin;


angle_step = 360 / args.num_frames;
profile_angles = 0:angle_step:(args.num_frames-1)*angle_step;

ylims = [min(input_patch(:)) max(input_patch(:))];

%Display patch
f1 = figure('color', args.fig_bg_color, ...
    'windowstyle', 'normal', 'units', 'normalized', 'position', [0 0 1 1]);
a1 = subplot(1,2,1); imgray(input_patch);
colormap(args.colormap);
axis off;

%Set up profiles plot
a2 = subplot(1,2,2); 
axis([1 2*args.profile_radius+1 ylims]);
axis off;
hold all;
plot([1 2*args.profile_radius], [ylims(1) ylims(1)], 'k', 'linewidth', 2);
plot([args.profile_radius+1 args.profile_radius+1], ylims, 'k', 'linewidth', 2);

%Get size of figure
set(f1, 'units', 'pixels');
fig_size = get(f1, 'position');
frame_rect = [50 50 fig_size(3:4)-100];

%Loop through each profile angle
for az = profile_angles
    
    %Get profile
    vx = args.profile_radius*cosd(az);
    vy = -args.profile_radius*sind(az);
    [px py p] = improfile(input_patch, xc + [-vx vx], yc + [-vy vy], 2*args.profile_radius+1, 'bilinear');
    
    %Update the plots
    if ~az
        plot(a1, px, py, [args.axes_color '--'], 'linewidth', 3,...
            'xdatasource', 'px', 'ydatasource', 'py');
        plot(a2, p, args.profile_color, 'linewidth', 3, 'ydatasource', 'p');
    else
        refreshdata(f1, 'caller');
    end
    
    %Get frames of figure and convert to gif format
    loop = true;
    while loop
        frame1 = getframe(f1, frame_rect);
        gif1 = frame2im(frame1);   
        [gif1a map] = rgb2ind(gif1, 2^16);
        loop = ~isa(gif1a, 'uint8');
    end
    
    %Save in gif
    if ~az        
        imwrite(gif1a, map, gif_name, 'gif', 'WriteMode', 'overwrite', 'Loop', Inf, 'DelayTime', args.delay_time);
    else        
        imwrite(gif1a, map, gif_name, 'gif', 'WriteMode', 'append', 'DelayTime', args.delay_time);
    end            
end