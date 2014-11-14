function [vessel_centre_i outer_edge_i inner_edge_i] = ...
                convert_edge_markup(outer_edge, inner_edge)
%CONVERT_EDGE_MARKUP *Insert a one line summary here*
%   [] = convert_edge_markup(contour_dir)
%
% Inputs:
%      contour_dir - *Insert description of input variable here*
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
% Created: 16-Apr-2013
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester

%Use spline function to get even spacing of points on both edges 
%(sampling much more densley on the inner edge) 
[inner_edge] = spline_contour(inner_edge, [], 0.1);
[outer_edge] = spline_contour(outer_edge, [], 1);

%Match each outer point to its nearest inner point (using only inner points
%approximately the same distance along the vessel and lying in a direction
%close to the normal projection to the outer edge at that point)
n_inner = size(inner_edge,1);
n_outer = size(outer_edge,1);
inner_ratios = (1:n_inner) / n_inner;
outer_ratios = (1:n_outer)/ n_outer;

%Compute normals to ecah outer edge point
outer_normals = compute_spline_normals(outer_edge);

outer_vec1 = inner_edge(1,:) - outer_edge(1,:);
outer_vec1 = outer_vec1 / sqrt(sum(outer_vec1.^2));

if (outer_vec1(:,1)*outer_normals(1,1) + outer_vec1(:,2)*outer_normals(1,2)) < 0
    outer_normals = -outer_normals;
end

%Set up containers for the matched inner points (and associated indices
%into the full set of inner edge pts)
matched_inner_edge = NaN(n_outer,2);
matched_inner_edge_idx = zeros(n_outer,1);
matched_outer_edge_idx = zeros(n_inner,1);
    
%Loop through each outer edge point
for i_pt = 1:n_outer

    %Compute distance and direction to each inner edge point
    outer_vecs = bsxfun(@minus, inner_edge, outer_edge(i_pt,:));
    dists = sum(outer_vecs.^2, 2);
    outer_vecs = bsxfun(@rdivide, outer_vecs, sqrt(dists));             

    %Select valid points
    local_inner_pts_idx = ...
        ( abs(inner_ratios - outer_ratios(i_pt)) < 0.1 )' & ...
        ( (outer_vecs(:,1)*outer_normals(i_pt,1) + outer_vecs(:,2).*outer_normals(i_pt,2)) > 0 );

    %If there's at least one valid point
    if any(local_inner_pts_idx)
        
        %Find the nearest and save its value and index in the appropiate
        %containers
        dists(~local_inner_pts_idx) = inf; 
        [~, nearest_idx] = min(dists);
        matched_inner_edge(i_pt, :) = inner_edge(nearest_idx,:);
        matched_inner_edge_idx(i_pt) = nearest_idx;
        matched_outer_edge_idx(nearest_idx) = ...
            matched_outer_edge_idx(nearest_idx) + 1;
    else
        0;
    end
end

%Keep only those poitns for which we have a matched pair
valid_pts = ~isnan(matched_inner_edge(:,1));
if ~all(valid_pts)
    outer_edge = outer_edge(valid_pts,:);
    matched_inner_edge = matched_inner_edge(valid_pts,:);
    matched_inner_edge_idx = matched_inner_edge_idx(valid_pts,:);
end  

if ~issorted(matched_inner_edge_idx)
    [matched_inner_edge_idx sort_idx] = sort(matched_inner_edge_idx);
    matched_inner_edge = matched_inner_edge(sort_idx,:);
end

%Take the centre of the two matched edges, and interpolate to evenly space
%the centre points
vessel_centre = (outer_edge + matched_inner_edge)/2;
[vessel_centre_i, dists, dists_i] = spline_contour(vessel_centre, [], 2, 'linear');  
num_pts = size(vessel_centre_i,1);

%Interpolate the outer points at the position of the evenly spaced centre
%points
outer_edge_i = interp1(dists, outer_edge, dists_i, 'linear');
inner_edge_i = zeros(num_pts,2);

%Now find an appropriately positioned inner edge for each centre/outer pair
%from the original set of inner edge points by projecting the vector from
%outer to centre and seeing where this crosses the inner edge
inner_edge_i(1,:) = inner_edge(1,:);
inner_edge_i(num_pts,:) = inner_edge(end,:);
    
for i_pt = 2:num_pts-1
    i1 = sum(dists < dists_i(i_pt));
    inner_idx1 = matched_inner_edge_idx(i1);
    inner_idx2 = matched_inner_edge_idx(i1+1);

    potential_pts = inner_edge(inner_idx1:inner_idx2,:);
    
    if isempty(potential_pts);
        0;
    end

    o_to_c = vessel_centre_i(i_pt,:) - outer_edge_i(i_pt,:);
    o_to_c = o_to_c / sqrt(sum(o_to_c.^2, 2));

    o_to_i = bsxfun(@minus, potential_pts, outer_edge_i(i_pt,:));
    o_to_id = sqrt(sum(o_to_i.^2, 2));
    o_to_i = bsxfun(@rdivide, o_to_i, o_to_id);

    dot_prods = o_to_i(:,1)*o_to_c(:,1) + o_to_i(:,2)*o_to_c(:,2);
    [~, nearest_idx] = max(dot_prods);

    inner_edge_i(i_pt,:) = inner_edge(inner_idx1 + nearest_idx - 1,:);
end

%Take the new centre as the average of the new inner and outer points
vessel_centre_i = (outer_edge_i + inner_edge_i)/2;
