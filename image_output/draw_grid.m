function [] = draw_grid(height,width,pos,step_size,color,axes_handles)
%DRAW_GRID *Insert a one line summary here*
%   [] = draw_grid(height,width,step_size,axes_handles,color)
%
% Inputs:
%      height- grid height
%
%      width- grid width
%
%      pos- [x y] coordinates of top left corner
%
%      step_size- size of grid squares
%
%      axes_handles- (optional - default current axes)
%
%      color- (optional - default 'r')
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 11-Apr-2008
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester 

if exist('axes_handles', 'var')
    axes(axes_handles);
end
if ~exist('step_size', 'var')
    step_size = 1;
end
if ~exist('color', 'var')
    color = 'r';
end

%Compute end points for grid lines
x_pts = pos(1):step_size:pos(1) + width;
y_pts = pos(2):step_size:pos(2) + height;

hold on;
%plot vertical grid lines (varying x)
for i = 1:length(x_pts)
    plot([x_pts(i) x_pts(i)], [y_pts(1) y_pts(end)], color);
end
%plot horizontal grid lines (varying y)
for i = 1:length(x_pts)
    plot([x_pts(1) x_pts(end)], [y_pts(i) y_pts(i)], color);
end