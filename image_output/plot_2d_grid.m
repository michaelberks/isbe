function [edges] = plot_2d_grid(xy, edges, edge_specs, edge_thickness, corner_specs, corner_size, plot_h)
%PLOT_CUBE *Insert a one line summary here*
%   [] = plot_cube(xyz, marker_type, edge_color, edge_thickness)
%
% Inputs:
%      xyz - *Insert description of input variable here*
%
%      edge_specs - *Insert description of input variable here*
%
%      edge_thickness - *Insert description of input variable here*
%
%      corner_specs - *Insert description of input variable here*
%
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 09-Mar-2017
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
if ~exist('xy', 'var') || isempty(xy)
    xy = [
        0 0  %Face 1
        1 0
        1 1
        0 1];
elseif size(xy,2) ~= 2
    error('Input xyz must be nx2 array, containing (x,y) coordinates of grid corners');
end

if ~exist('edge_specs', 'var') || isempty(edge_specs)
    edge_specs = '-';
end
if ~exist('edge_thickness', 'var') || isempty(edge_thickness)
    edge_thickness = 2;
end
if ~exist('corner_specs', 'var') || isempty(corner_specs)
    corner_specs = 'b*';
end
if ~exist('corner_size', 'var') || isempty(corner_size)
    corner_size = 10;
end

%Each of x,y,z lines is mx2 array defining the start and end of each edge
%in the respective dimension

if ~exist('edges', 'var') || isempty(edges)
    
    edges = zeros(0,2);
    num_pts = size(xy, 1);
    for i_pt = 1:num_pts-1
        for j_pt = i_pt+1:num_pts
            %Check if adjacent
            if sum(abs(xy(i_pt,:)-xy(j_pt,:))) < 1.1
                edges = [edges; i_pt j_pt]; %#ok
            end
        end
    end
    if strcmpi(edge_specs, 'no_plot')
        return;
    end
end

if ~exist('plot_h', 'var') || isempty(plot_h)
    plot_h = gca;
end

if numel(plot_h) == 1
    plot(plot_h, ...
        [xy(edges(:,1),1) xy(edges(:,2),1)]',...
        [xy(edges(:,1),2) xy(edges(:,2),2)]',...
        edge_specs, 'LineWidth', edge_thickness); hold on;
    plot(plot_h, xy(:,1), xy(:,2), corner_specs, 'MarkerSize', corner_size);

end


