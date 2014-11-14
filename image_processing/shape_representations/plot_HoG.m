function [axes_h] = plot_HoG(image_in, hog, cell_sz, angle_wrap, axes_h)
%COMPUTE_HOG *Insert a one line summary here*
%   [] = compute_HoG(varargin)
%
% COMPUTE_HOG uses the U_PACKARGS interface function
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
% Created: 30-Sep-2013
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 


%Make cell x,y coords
[n_rows n_cols] = size(image_in);
cl_x = (cell_sz(2) / 2):cell_sz(2):n_cols;
cl_y = (cell_sz(1) / 2):cell_sz(1):n_rows;

%Check we have the correct number of cells
[num_cells num_ori_bins] = size(hog);
n_cl_rows = length(cl_y);
n_cl_cols = length(cl_x);
if (n_cl_rows*n_cl_cols) ~= num_cells
    error('Incorrect number of cells in HoG representation');
end

%Make theta vectors for each orientation bin
ori_theta = fliplr(linspace(0, 2*pi, num_ori_bins+1));
if angle_wrap
    ori_theta = ori_theta / 2;
end
ori_theta = (ori_theta(1:end-1)+ori_theta(2:end)) / 2; 
ori_x = cos(ori_theta);
ori_y = -sin(ori_theta);

%Normalise the HoG scores and get the smaller of the cell height/width so
%we make sure no line overlaps the next cell in the plot
max_hog = max(hog(:));
hog_n = hog / max_hog;
min_cell_dim = min(cell_sz);

%Draw image in
if ~exist('axes_h', 'var') || isempty(axes_h)
    figure;
    axes_h = gca;
else
    axes(axes_h); %#ok
end
imgray(image_in);
cell_num = 1;

%Loop through the cells
for i_col = 1:n_cl_cols
    cl_x_i = cl_x(i_col);
    
    for i_row = 1:n_cl_rows
        
        cl_y_i = cl_y(i_row);
        
        for i_ori = 1:num_ori_bins
            
            %Plot each HoG score as line, oriented to the bin direction,s
            %caled by the bin count
            length_i = hog_n(cell_num, i_ori);
            hog_x_i = cl_x_i + length_i*ori_x(i_ori)*[-1 1]*min_cell_dim/2;
            hog_y_i = cl_y_i + length_i*ori_y(i_ori)*[-1 1]*min_cell_dim/2;
            
            plot(axes_h, hog_x_i, hog_y_i, 'r');
        end
        cell_num = cell_num+1;
    end
end
            
            
            

    
    
    




