function [x_pts, y_pts, selected_markers, marker_ui_xy] = ...
    marker_auto_detect(mammo, max_pairs, marker_template, large_mammo, left_breast, reduction_factor, debug_mode)
%MARKERDETECT detect breast thickness markes in mammogram
%
% Arguments:
%
% 'mammo'
%   - Mammogram containing breast thickness markers
%
% 'max_pairs'
%   - Maximum number of pairs user will be asked to select
%
% 'debug_mode (0)'
%   - turn image output on/off
%
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

if nargin < 3
    debug_mode = 0;
end
rows = size(mammo,1);
show_rows = 800;
missed_portion = rows - 2*show_rows;

%Apply the auto detection algorithm
[markers region_lims] = ...
    find_markers(mammo ,marker_template, large_mammo, left_breast, reduction_factor, debug_mode, debug_mode);
region_centres = [mean(region_lims(:,[1 2]),2) mean(region_lims(:,[3 4]),2)];

num_markers = size(markers,1);
if num_markers ~= 2*max_pairs
    error('Incorrect number of markers returned from FIND_MARKERS');
end

show_lims = region_lims;
show_lims(max_pairs+1:end,[3 4]) = show_lims(max_pairs+1:end,[3 4]) - missed_portion;
show_centres = region_centres;
show_centres(max_pairs+1:end,2) = show_centres(max_pairs+1:end,2) - missed_portion;
show_markers = markers(:,[1 2]);
show_markers(max_pairs+1:end,2) = show_markers(max_pairs+1:end,2) - missed_portion;

%Display image, only show the first and last 600 rows so the markers are
%more visible

marker_fig = figure(...
    'Name', 'Review marker points',...
    'Units', 'normalized',...
    'OuterPosition', [0 0 1 1]);
imgray(mammo([1:show_rows rows-show_rows+1:rows],:));

for i_m = 1:2*max_pairs %Should also be #rows in markers   
        plot(show_lims(i_m,[1 2 2 1 1]), show_lims(i_m,[3 3 4 4 3]));
        if markers(i_m,1)
            plot(show_markers(i_m,1), show_markers(i_m,2), 'go', 'MarkerSize',10);
            plot(show_markers(i_m,1), show_markers(i_m,2), 'gx', 'MarkerSize',10);
        else
            text(show_centres(i_m,1), show_centres(i_m,2), 'Missing', 'color', 'r');
        end
end

title({['Have all ' num2str(max_pairs) ' pairs been correctly found or identified as missing?'];
    'If ''yes'' press enter, otherwise click on each marker region that needs correcting, then press enter'});
set(gca, 'xticklabel', [], 'yticklabel', []);

[xi,yi,P] = impixel; %#ok

if isempty(xi)
    marker_ui_xy = [];
else
    marker_ui_xy = [xi(:) yi(:)];

    for i_m = 1:length(xi);
        d = (show_centres(:,1)-xi(i_m)).^2 + (show_centres(:,2)-yi(i_m)).^2;
        [~,r] = min(d);
        region_c = region_lims(r,1):region_lims(r,2);
        region_r = region_lims(r,3):region_lims(r,4);

        %Extract region from image and mask
        mammo_i = mammo(region_r, region_c);

        [x y] = correct_marker(mammo_i);
        if isempty(x)
            markers(r,1) = 0;
            markers(r,2) = 0;
        else
            markers(r,1) = x + region_lims(r,1) - 1;
            markers(r,2) = y + region_lims(r,3) - 1;
        end
    end
end

%Now sort out which pairs we have
% x_pts/y_pts will hold the coords of the centre of each marker
x_pts = zeros(2, max_pairs);
y_pts = zeros(2, max_pairs);

for i_p = 1:max_pairs
    x_pts(1,i_p) = markers(i_p,1);
    x_pts(2,i_p) = markers(i_p+max_pairs,1);
    
    y_pts(1,i_p) = markers(i_p,2);
    y_pts(2,i_p) = markers(i_p+max_pairs,2);
end

selected_markers = find(all(x_pts));
x_pts = x_pts(:,selected_markers);
y_pts = y_pts(:,selected_markers);

close(marker_fig);