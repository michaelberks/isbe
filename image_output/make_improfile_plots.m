function [] = make_improfile_plots(input_patch, im_name, xc, yc, varargin)
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
    'profile_angles', [-75 -45 -15 15 45 75],...
    'fig_bg_color', 'white',...
    'image_height', 400,...
    'ylims', [],...
    'colormap', gray(256));
clear varargin;

if isempty(args.ylims)
    ylims = [min(input_patch(:)) max(input_patch(:))];
else
    ylims = args.ylims;
end

[rows cols] = size(input_patch);
im_h = args.image_height;
im_w = im_h*cols/rows;

%Display patch
f1 = figure('color', args.fig_bg_color, ...
    'windowstyle', 'normal', 'units', 'pixels', 'position', [50 50 im_w im_h]);
a1 = axes('units', 'normalized', 'position', [0 0 1 1]);
imgray(input_patch);
colormap(args.colormap);
axis off;

%Set up profiles plot
f2 = figure('color', args.fig_bg_color, ...
    'windowstyle', 'normal', 'units', 'pixels', 'position', [450 50 400 400]);
a2 = axes('units', 'normalized', 'position', [0 0 1 1]);
axis([1 2*args.profile_radius+1 ylims]);
axis off;
hold all;
plot([1 2*args.profile_radius], [ylims(1) ylims(1)], 'k', 'linewidth', 2);
plot([args.profile_radius+1 args.profile_radius+1], ylims, 'k', 'linewidth', 2);


profile_colors = lines(length(args.profile_angles));

%Loop through each profile angle
for i_az = 1:length(args.profile_angles);
    
    %Get profile
    az = args.profile_angles(i_az);
    vx = args.profile_radius*cosd(az);
    vy = -args.profile_radius*sind(az);
    [px py line_profile] = improfile(input_patch, xc + [-vx vx], yc + [-vy vy], 2*args.profile_radius+1, 'bilinear');
    
    plot(a1, px, py, '--', 'color', profile_colors(i_az,:), 'linewidth', 3);
    %Update the plots
    if i_az == 1        
        p1 = plot(a2, line_profile, 'color', profile_colors(1,:), 'linewidth', 3, 'ydatasource', 'line_profile');
    else
        set(p1, 'color', profile_colors(i_az,:)); 
        refreshdata(f2, 'caller');        
    end
    
    figure(f2);
    exportfig([im_name '_' zerostr(i_az,2) '.png']);
    delete([im_name '_' zerostr(i_az,2) '.fig']);
    delete([im_name '_' zerostr(i_az,2) '.pdf']);
    delete([im_name '_' zerostr(i_az,2) '.eps']);         
end
% figure(f1);
% exportfig([im_name '.png']);
% delete([im_name '.fig']);
% delete([im_name '.pdf']);
% delete([im_name '.eps']);

print('-dtiff', '-noui', '-painters', f1, '-r300', [im_name '.tif']);