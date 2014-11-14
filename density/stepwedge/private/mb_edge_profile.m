function [breast_thickness taper_region_idx] = mb_edge_profile(x_b, breast_border, breast_air, left_breast, resize_factor, mlo, debug_mode)
%MB_EDGE_PROFILE compute thickness of breats as it tapers at breast edge
%
% Arguments:
%
% 'x_b' - map of breast thickness
%
% 'breast_border' - coordinates of breast edge
%
% 'breast_air' - indices to border coordinate son breast/air boundary
%
% 'left_breast' - Flag: 1 if left, 0 if right
%
% 'Resize factor' - Image resolution...
%
% 'mlo' - Flag: 1 if MLO, 0 if CC
%
% Outputs: none
%
% 'x_pts, y_pts'
%   - x,y coordinates of select markers
%
% 'r_pts'
%   - radius of each marker
%
% 'selected_markers'
%   - Binary vector of length max_pairs, indicating whether marker pair was
%   present (1) or missing (0)
%
% 'errorcheck'
%  - return 0 unless user has selected no marker pairs are visible in the
%  image
%
% See also: STEPWEDGE, FIND_MARKER
%
% Created: 09-May-2006
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 5125 
% Copyright: (C) University of Manchester

if nargin < 7
    debug_mode = 0;
end

%Get dimensions of mammogram
dims = size(x_b);

%If right mammogram flip thickness map so chestwall is on left
if ~left_breast
    x_b=flipdim(x_b,2);
    
    % flip border points too
    breast_border(:,1) = dims(2) - breast_border(:,1) + 1;
    
end

if debug_mode
    figure;
    imagesc(x_b); axis image; colormap(gray(256)); hold on;
    plot(breast_border(:,1), breast_border(:,2), 'c', 'LineWidth', 2);
end
 
%Get length of edge boundary
numedgepoints = length(breast_air);

%Pre-allocate space for profiles and inner edge
profile_length = 50; %MB have this as an argument to set
x_profiles = zeros(numedgepoints, profile_length);
y_profiles = zeros(numedgepoints, profile_length);
thickness_profiles = zeros(numedgepoints, profile_length);
inner_xy = zeros(numedgepoints, 2);
    
%Move through edge-points calculating thicknesses along normal profiles
%from breast edge inwards, assuming a semi-circular profile;
for ii = 1:numedgepoints
    
    centrex = breast_border(breast_air(ii),1);
    centrey = breast_border(breast_air(ii),2);
    
    si = max(1, ii-5);
    ei = min(numedgepoints, ii+5);
    
    %Find unit vector between adjacent edge points of centre point
    tan_vec = [breast_border(breast_air(ei),1) - breast_border(breast_air(si),1)...
               breast_border(breast_air(ei),2) - breast_border(breast_air(si),2)];
    tan_vec = tan_vec / (sqrt(sum(tan_vec.^2)));
    
    %Convert tangent vector to normal vector
    normvec = [-tan_vec(2) tan_vec(1)];

    % get thickness at the centre point, divide by 2 to get radius - need to
    % correct to pixels not mm!!!
    centrethick = double(x_b(round(centrey),round(centrex)))*1e-3*resize_factor/(44e-6);
    centrethick = centrethick / 2.;
    
    % calculate point on the normal which is centrethick away from the edge
    inner_xy(ii, :) = [centrex centrey] + centrethick*normvec;
    
    if debug_mode
        plot([inner_xy(ii, 1) centrex], [inner_xy(ii, 2) centrey], 'g');
        plot(centrex, centrey, 'rx');
    end
    
    % Now calculate x-y profiles moving from the edge to a point
    % centrethick away
    pts = linspace(0, centrethick, profile_length);
    x_profiles(ii, :) = centrex + pts*normvec(1);
    y_profiles(ii, :) = centrey + pts*normvec(2);
    
    %Finally calculate the thickness based on a elliptical profile at
    %each point - the ellipse has a horizontal radius of centrethick, and a
    %vertical radius equal to the thickness of the breast at the inner
    %point on the profile
    r = max(1, round(inner_xy(ii,2)));
    c = max(1, round(inner_xy(ii,1)));
    
    %Note there is a factor of 2 in innerthick that comes out in the wash
    %when we apply the profile - so don't divide by 2 like centrethick
    innerthick = x_b(r,c);

    %Reverse points so max of profile is at inner point, 0 at centre point
    pts = pts(end:-1:1);
    thickness_profiles(ii, :) = ...
        innerthick * sqrt(1 - (pts.^2 / centrethick^2));
    
end

%Extend the inner border to the end of the image (depending on view type)
if mlo
    %if MLO at start fit x at top of mam, at end fit y to chestwall
    p_start = polyfit(inner_xy(1:5,2),inner_xy(1:5,1),1);
    p_end = polyfit(inner_xy(end-4:end,1),inner_xy(end-4:end,2),1);
    inner_xy = [sum(p_start) 1; inner_xy; 1 sum(p_end)];
    
else
    %if CC at start and end fit y to chestwall
    p_start = polyfit(inner_xy(1:5,1), inner_xy(1:5,2),1);
    p_end = polyfit(inner_xy(end-4:end,1), inner_xy(end-4:end,2),1);
    inner_xy = [1 sum(p_start); inner_xy; 1 sum(p_end)];
end

% Compute thickness profiles for the new start and end points
%start point
centrex = breast_border(breast_air(1),1);
centrey = breast_border(breast_air(1),2);

centrethick = double(x_b(round(centrey),round(centrex)))*1e-3*resize_factor/(44e-6);
centrethick = centrethick / 2.;

dirvec = inner_xy(1,:) - [centrex centrey];
dirvec = dirvec / sqrt(sum(dirvec.^2));
pts = linspace(0, centrethick, profile_length);
x_profiles = [centrex + pts*dirvec(1); x_profiles];
y_profiles = [centrey + pts*dirvec(2); y_profiles];
r = max(1, round(inner_xy(1,2)));
c = max(1, round(inner_xy(1,1)));
innerthick = x_b(r,c);
pts = pts(end:-1:1);
thickness_profiles = [innerthick * sqrt(1 - (pts.^2 / centrethick^2)); thickness_profiles];

%end point
centrex = breast_border(breast_air(end),1);
centrey = breast_border(breast_air(end),2);

centrethick = double(x_b(round(centrey),round(centrex)))*1e-3*resize_factor/(44e-6);
centrethick = centrethick / 2.;

dirvec = inner_xy(end,:) - [centrex centrey];
dirvec = dirvec / sqrt(sum(dirvec.^2));
pts = linspace(0, centrethick, profile_length);
x_profiles = [x_profiles; centrex + pts*dirvec(1)];
y_profiles = [y_profiles; centrey + pts*dirvec(2)];
r = max(1, round(inner_xy(end,2)));
c = max(1, round(inner_xy(end,1)));
innerthick = x_b(r,c);
pts = pts(end:-1:1);
thickness_profiles = [thickness_profiles;
    innerthick * sqrt(1 - (pts.^2 / centrethick^2))];


%Reverse the outer border to create a polygon bounding the area between the
%inner and outer border
taper_border = [inner_xy; flipud(breast_border(breast_air,:))];

%Create mask of inner area
breast_mask = poly2mask(breast_border(:,1), breast_border(:,2), dims(1), dims(2));
taper_mask = poly2mask(taper_border(:,1), taper_border(:,2), dims(1), dims(2));
taper_mask = taper_mask & breast_mask;



%Get x,y-coordinates of region to fill (i.e. difference between inner and
%outer area
[taper_region_y taper_region_x] = find(taper_mask);
taper_region_idx = sub2ind(dims, taper_region_y, taper_region_x);

%Use griddata to fill in the missing data
thicknesses = griddata(x_profiles(:), y_profiles(:), thickness_profiles(:), taper_region_x, taper_region_y, 'linear');
thicknesses(isnan(thicknesses)) = 0;

%Now put all this together into final thickness map
breast_thickness = x_b;
breast_thickness(~breast_mask) = 0;
breast_thickness(taper_region_idx) = thicknesses;

%flip everything back
if ~left_breast
    breast_thickness = flipdim(breast_thickness, 2);
    taper_region_x = dims(2) + 1 - taper_region_x;
    taper_region_idx = sub2ind(dims, taper_region_y, taper_region_x);
end

%MB only show in debug_mode
if debug_mode
    figure; imagesc(breast_thickness); axis image; colormap(gray(256));
    figure; imagesc(taper_mask); axis image;  colormap(gray(256));
end





