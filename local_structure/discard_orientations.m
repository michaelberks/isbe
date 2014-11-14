function [discard_mask] = discard_orientations(ori_map, varargin)
%DISCARD_ORIENTATIONS *Insert a one line summary here*
%   [] = discard_orientations(varargin)
%
% DISCARD_ORIENTATIONS uses the U_PACKARGS interface function
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
% Created: 23-Mar-2011
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

% Unpack the arguments:
args = u_packargs(varargin, '0', ...
    'angle_tol', pi/24,...
    'mask', [],...
    'edge_size', 50);
clear varargin;
edge_size = args.edge_size;
angle_tol = args.angle_tol;
mask = args.mask;
clear args;

[row col] = size(ori_map);

%Create masks of angles to discard at top, bottom, left and right
left_mask = [true(row, edge_size) false(row, col-edge_size)] &...
    (mod(ori_map,pi) > (pi/2 - angle_tol)) &...
    (mod(ori_map,pi) < (pi/2 + angle_tol));

right_mask = [false(row, col-edge_size) true(row,edge_size)] &...
    (mod(ori_map,pi) > (pi/2 - angle_tol)) &...
    (mod(ori_map,pi) < (pi/2 + angle_tol));

top_mask = [true(edge_size, col); false(row-edge_size, col)] &...
    ((mod(ori_map,pi) > (pi-angle_tol)) | (mod(ori_map,pi) < angle_tol));


bottom_mask = [false(row-edge_size, col); true(edge_size, col)] &...
    ((mod(ori_map,pi) > (pi-angle_tol)) | (mod(ori_map,pi) < angle_tol));

discard_mask = left_mask | right_mask | top_mask | bottom_mask;

if ~isempty(mask)
    
    %trace the boundary of the map and convert to xy form
    [st_r st_c] = find(mask, 1);
    xy = bwtraceboundary(mask, [st_r st_c], 'N', 4, inf, 'clockwise');
    xy = fliplr(xy(1:end-1,:)); %Last point is repeat of first
    num_pts = size(xy,1);
    
    %wrap the start and end of the boundary
    xy_wrap = [xy(end-9:end,:); xy; xy(1:10,:)];
    
    %Get directions at each point on the boundary
    uv = zeros(num_pts,2);
    for ii = 1:num_pts
        i1 = ii;
        i2 = ii+20;
        uv(ii,:) = mean(diff(xy_wrap(i1:i2,:)));
        uv(ii,:) = uv(ii,:) / sqrt(sum(uv(ii,:).^2));
    end
    
    %Now make a mask of the edge region
    edge_mask = mask & ~imerode(mask, strel('disk', edge_size));
    
    %Get xy coordinates of the inner mask
    [yi xi] = find(edge_mask);

    %Use griddata to interpolate the orientations from the boundary edge to
    %the whole of the edge mask
    uvi = zeros(length(yi),2);
    uvi(:,1) = griddata(xy(:,1),xy(:,2),uv(:,1),xi,yi, 'nearest');
    uvi(:,2) = griddata(xy(:,1),xy(:,2),uv(:,2),xi,yi, 'nearest');

    %Copy these orientations into image form
    edge_ori = zeros(row, col);
    edge_ori(edge_mask) = atan2(-uvi(:,2), uvi(:,1));
    
    %Workout orientations that align with the edge mask
    edge_mask = edge_mask & abs(mb_mod(ori_map - edge_ori,pi)) < (2*angle_tol);
    
    %Combine with the discard mask
    discard_mask = discard_mask | edge_mask;

end


    
    
        
    