function [x_pts, y_pts, r_pts, selected_markers, marker_ui_xy, errorcheck] = markerdetect(mammo, max_pairs, left_breast, debug_mode)
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

%Display image, only show the first and last 600 rows so the markers are
%more visible
show_rows = 600;
marker_fig = figure('Name',['Select marker points on this image - choose at most ', num2str(max_pairs), ' pairs']);
imagesc(mammo([1:show_rows rows-show_rows+1:rows],:)); axis image; colormap(gray(256));
missed_portion = rows - 2*show_rows;

if max_pairs == 3
    title('Choose 3 pairs, ignoring the pair furthest to nipple side');
elseif max_pairs == 4
    title('Choose 4 pairs');
end
if left_breast
    set(gca, 'yaxislocation', 'right');
end
xlabel('Press enter when finished');
ylabel('Nipple side');
set(gca, 'xticklabel', [], 'yticklabel', []);

[xi,yi,P] = impixel; %#ok
marker_ui_xy = [xi(:) yi(:)];

% find how many markers were selected
numpoints = length(xi);
errorcheck = 0;

%if 0 or 1 points skip this image
if numpoints < 2
    answer = questdlg({'You have selected 1 or fewer markers.';...
        'If this is because no markers are present, click Yes to continue to the next mammogram.';...
        'Otherwise, click No to select more points,'}, '0 or 1 points selected','Yes','No','Yes');
    
    close(marker_fig);
    
    if strcmpi(answer, 'yes');
        errorcheck = 1;
        x_pts = 0;
        y_pts = 0;
        r_pts = 0;
        selected_markers = 0;
        return;
    else
        [x_pts,y_pts,r_pts,selected_markers,marker_ui_xy,errorcheck] = markerdetect(mammo, max_pairs, debug_mode);
        return
    end

%if odd number of points or too many points, try again
elseif rem(numpoints, 2)
    close(marker_fig);
    warndlg('Odd number of points selected, please start again selecting only complete pairs', 'Wrong number of points selected', 'modal');
    [x_pts,y_pts,r_pts,selected_markers,marker_ui_xy,errorcheck] = markerdetect(mammo, max_pairs, debug_mode);
    return
elseif max_pairs == 3 && numpoints > 6
    close(marker_fig);
    warndlg({'Ignoring the pair closest to the nipple means there should be at most three pairs'; 'Please start again'}, 'Too many points selected', 'modal');
    [x_pts,y_pts,r_pts,selected_markers,marker_ui_xy,errorcheck] = markerdetect(mammo, max_pairs, debug_mode);
    return
elseif max_pairs == 4 && numpoints > 8
    close(marker_fig);
    warndlg({'There should be at most four pairs'; 'Please start again'}, 'Too many points selected', 'modal');
    [x_pts,y_pts,r_pts,selected_markers,marker_ui_xy,errorcheck] = markerdetect(mammo, max_pairs, debug_mode);
    return
end
%Otherwise keep going

%Work out number of pairs
num_pairs = numpoints / 2;

%Don't assume the points have been selected in the correct order, sort into
%pairs now
xy_pts = [xi(:) yi(:)];

%first sort by yi, to separate upper markers from lower markers
xy_pts = sortrows(xy_pts, 2);

%then sort each set of markers by xi
xy_pts(1:num_pairs,:) = sortrows(xy_pts(1:num_pairs,:), 1);
xy_pts(num_pairs+1:numpoints,:) = sortrows(xy_pts(num_pairs+1:numpoints,:), 1);

% x_pts/y_pts will hold the coords of the centre of each marker
x_pts = zeros(num_pairs, 2);
y_pts = zeros(num_pairs, 2);
r_pts = zeros(num_pairs, 2);

% loop over the marker pairs finding each marker
for ii = 1:numpoints
    %Detect marker automatically
    [x_pts(ii), y_pts(ii), r_pts(ii)] =...
        find_marker(mammo([1:show_rows rows-show_rows+1:rows],:), xy_pts(ii,1), xy_pts(ii,2), 0, marker_fig);
    
end

%Transpose points so pairs are in columns
x_pts = x_pts';
y_pts = y_pts';
r_pts = r_pts';

%we also need to account for the missed y-portion of the image in the bottom
%markers
y_pts(2,:) = y_pts(2,:) + missed_portion;

%If less than max_pairs selected need to find out which pairs have been
%selected
selected_markers = zeros(max_pairs,1);

if num_pairs < max_pairs % in this case need to sort out which marker pairs we have
                    % so that correct calibration factors are used
    for ii = 1:max_pairs
        answer = questdlg(['Is pair ', num2str(ii), ' present?'],'Marker Points','Yes','No','Yes');
        selected_markers(ii) = strcmpi(answer, 'yes');
    end
else
    selected_markers(:) = 1;
end
selected_markers = find(selected_markers);

close(marker_fig);