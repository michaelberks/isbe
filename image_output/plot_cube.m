function [] = plot_cube(xyz, edge_specs, edge_thickness, corner_specs, corner_size, plot_h)
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
if isempty(xyz)
    xyz = [
        0 0 0 %Face 1
        1 0 0
        1 1 0
        0 1 0
        0 0 1 %Face 2
        1 0 1
        1 1 1
        0 1 1];
elseif size(xyz,1) ~= 8 || size(xyz,2) ~= 3
    error('Input xyz must be 8x3 array, containing coordinates of cube corners');
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

if ~exist('plot_h', 'var') || isempty(plot_h)
    plot_h = gca;
end

%Each of x,y,z lines is 12x2 array defining the start and end of each edge
%in the respective dimension
x_edges = [
    %Top face
    xyz(1,1) xyz(2,1)
    xyz(2,1) xyz(3,1)
    xyz(3,1) xyz(4,1)
    xyz(1,1) xyz(4,1)
    %Bottom face
    xyz(5,1) xyz(6,1)
    xyz(6,1) xyz(7,1)
    xyz(7,1) xyz(8,1)
    xyz(5,1) xyz(8,1)
    %Edges joining two faces
    xyz(1,1) xyz(5,1)
    xyz(2,1) xyz(6,1)
    xyz(3,1) xyz(7,1)
    xyz(4,1) xyz(8,1)]';
y_edges =[
    %Top face
    xyz(1,2) xyz(2,2)
    xyz(2,2) xyz(3,2)
    xyz(3,2) xyz(4,2)
    xyz(1,2) xyz(4,2)
    %Bottom face
    xyz(5,2) xyz(6,2)
    xyz(6,2) xyz(7,2)
    xyz(7,2) xyz(8,2)
    xyz(5,2) xyz(8,2)
    %Edges joining two faces
    xyz(1,2) xyz(5,2)
    xyz(2,2) xyz(6,2)
    xyz(3,2) xyz(7,2)
    xyz(4,2) xyz(8,2)]';
z_edges = [
    %Top face
    xyz(1,3) xyz(2,3)
    xyz(2,3) xyz(3,3)
    xyz(3,3) xyz(4,3)
    xyz(1,3) xyz(4,3)
    %Bottom face
    xyz(5,3) xyz(6,3)
    xyz(6,3) xyz(7,3)
    xyz(7,3) xyz(8,3)
    xyz(5,3) xyz(8,3)
    %Edges joining two faces
    xyz(1,3) xyz(5,3)
    xyz(2,3) xyz(6,3)
    xyz(3,3) xyz(7,3)
    xyz(4,3) xyz(8,3)]';

plot3(plot_h, x_edges, y_edges, z_edges, edge_specs, 'LineWidth', edge_thickness); hold on;
plot3(plot_h, xyz(:,1), xyz(:,2), xyz(:,3), corner_specs, 'MarkerSize', corner_size);


