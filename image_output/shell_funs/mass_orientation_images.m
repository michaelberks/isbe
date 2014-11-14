function [] = mass_orientation_images(ii, ff, varargin)
%TEST_ORIENTATION_IMAGES *Insert a one line summary here*
%   [] = test_orientation_images(ii, ff)
%
% Inputs:
%      ii - *Insert description of input variable here*
%
%      ff - *Insert description of input variable here*
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
% Created: 11-Jan-2011
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

args = u_packargs(varargin,... % the user's input
    '0', ... % non-strict mode
    'roi_list', [],...
    'roi_dir', [asymmetryroot 'data/mammograms/2004_screening_processed/mass_roi/'],...
    'line_dir', [asymmetryroot 'data/line_maps/g2d/2004_screening_processed/mass_roi/'],...
    'orientation_dir', [asymmetryroot 'data/orientation_maps/g2d/2004_screening_processed/mass_roi/'],...
    'meta_dir', [asymmetryroot 'data/mammograms/2004_screening_processed/mass_roi/meta/'],...
    'spacing', 2,...
    'num_angles', 36,...
    'quiver_length', [],...
    'line_thresh', 0);
clear varargin;

ang_res = 180 / args.num_angles;
arrow_colors = hsv(args.num_angles);

figure(ff);

warning('off', 'load_uint8:missing_variables');

%If not supplied, get list of mammogram rois
if isempty(args.roi_list)
    roi_list = dir([args.roi_dir '*.mat']);
end

%Get mammo names
[roi_names] = get_mammo_info(roi_list);
num_images = length(roi_list);

%Get list of line rois
[line_names missing_idx] =...
    match_mammo_names(args.line_dir, roi_names);

%Get list of line rois
[orientation_names orientation_idx] =...
    match_mammo_names(args.orientation_dir, roi_names);
missing_idx = union(missing_idx, orientation_idx);

%Get list of mass meta info
[meta_names meta_idx] =...
    match_mammo_names(args.meta_dir, roi_names);
missing_idx = union(missing_idx, meta_idx);

%Now work out valid images
if isempty(missing_idx)
    valid_idx = 1:num_images;
else
    warning('mass_feature:incomplete','There was incomplete data for this set of mammograms');
    valid_idx = setdiff(1:num_images, missing_idx)';
end

ii = valid_idx(ii);

roi = load_uint8([args.roi_dir roi_list(ii).name]);
line_map = load_uint8([args.line_dir line_names{ii}]);
orientation_map = load_uint8([args.orientation_dir orientation_names{ii}]);
if ~isreal(orientation_map)
    orientation_map = 180*angle(orientation_map) / pi;
end
orientation_map = mod(orientation_map,180);
[r c] = size(orientation_map);
if any(size(line_map) ~= [r c])
    line_map = imresize(line_map, [r c], 'bilinear');
end
if any(size(roi) ~= [r c])
    roi = imresize(roi, [r c], 'bilinear');
end
mass_xy = u_load([args.meta_dir meta_names{ii}]);
mass_x = mass_xy(:,1)*c;
mass_y = mass_xy(:,2)*r;

a1 = subplot(1,2,1);
imagesc(roi); colormap(gray(256)); axis image; hold on;
plot(mass_x, mass_y, 'r:');
title(['Mass ROI ' roi_list(ii).name]);

a2 = subplot(1,2,2);
imagesc(line_map); colormap(gray(256)); axis image; hold on;
title('Line map with orientation quivers');


mask = line_map > args.line_thresh;

spacing_mask = false(r, c);
spacing_mask(1:args.spacing:r, 1:args.spacing:c) = true;
mask = mask & spacing_mask;

if isempty(args.quiver_length)
    args.quiver_length = r/100;
end

for jj = 0:args.num_angles
    theta = (jj - 0.5)*ang_res;

    %Get mask of pixels that have orientation within theta range
    angle_mask = mask & ...
         (orientation_map > theta - 0.5*ang_res) &...
         (orientation_map <= theta + 0.5*ang_res);

    [y x] = find(angle_mask);
    u = cos(pi*orientation_map(angle_mask)/180);
    v = -sin(pi*orientation_map(angle_mask)/180);

    quiver(a2, x, y, args.quiver_length*u, args.quiver_length*v, 0,...
        'color', arrow_colors(mod(jj, args.num_angles)+1,:));
end
plot(mass_x, mass_y, 'r');

linkaxes([a1 a2]);
zoom on;

hManager = uigetmodemanager(ff);
set(hManager.WindowListenerHandles,'Enable','off');
