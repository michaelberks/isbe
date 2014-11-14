function [] = compute_rad_maps_batch(data_type, varargin)
%COMPUTE_RAD_MAPS_BATCH *Insert a one line summary here*
%   [] = compute_rad_maps_batch()
%
% Inputs:
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 12-Jul-2010
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester

if ispc
    %task and job identifiers
    task_id = 1;
    num_jobs = 20;
    view = [];
    
    %location of data maps
    line_dir = 'line_maps';
    ori_dir = 'orientation_maps';
    relevance_dir = [];
    radial_dir = 'radial_maps';

    %arguments controlling computation of radial maps
    angular_resolution = 36;
    angular_bands = 1;
    distance_ranges = [128 256 inf];
    resize = 0.5;
    
else
%Set defaults for unix system - environment variables listed in .bashrc (and utils/mb_bash.m)
    
    %task and job identifiers
    [z task_id] = unix('echo $SGE_TASK_ID'); task_id = str2num(task_id);
    [z num_jobs] = unix('echo $NUM_JOBS'); num_jobs = str2num(num_jobs);
    [z view] = unix('echo $VIEW'); view(end) = [];
    
    %location of data maps
    [z line_dir] = unix('echo $LINE_DIR'); line_dir(end) = [];
    [z ori_dir] = unix('echo $ORI_DIR'); ori_dir(end) = [];
    [z relevance_dir] = unix('echo $RELEVANCE_DIR'); relevance_dir(end) = [];
    [z radial_dir] = unix('echo $RADIAL_DIR'); radial_dir(end) = [];
    
    %arguments controlling computation of radial maps
    [z angular_resolution] = unix('echo $ANGULAR_RESOLUTION'); angular_resolution = str2num(angular_resolution);
    [z angular_bands] = unix('echo $ANGULAR_BANDS'); angular_bands = str2num(angular_bands);
    [z distance_ranges] = unix('echo $DISTANCE_RANGE'); distance_ranges = str2num(distance_ranges);
    [z resize] = unix('echo $RESIZE'); resize = str2num(resize);

    clear z;
end

args = u_packargs(varargin,... % the user's input
    '0', ... % non-strict mode
    'task_id', task_id, ...
    'num_jobs', num_jobs, ...
    'view', view,...
    'line_dir', line_dir, ...
    'ori_dir', ori_dir, ...
    'relevance_dir', relevance_dir, ...
    'radial_dir', radial_dir, ...
    'angular_resolution', angular_resolution, ...
    'angular_bands', angular_bands, ...
    'distance_ranges', distance_ranges, ...
    'resize', resize);

angular_resolution = args.angular_resolution;
angular_bands = args.angular_bands;
distance_ranges = args.distance_ranges;
resize = args.resize;
    
%Get list of line maps
line_dir = [asymmetryroot 'data/' args.line_dir '/' data_type '/'];
line_list = dir([line_dir '*' args.view '*.mat']);
mam_names = get_mammo_info(line_list);

%Get list of orientation maps
ori_dir = [asymmetryroot 'data/' args.ori_dir '/' data_type '/'];

%Get list of relevance classification maps (if included)
if ~isempty(args.relevance_dir)
    relevance_dir = [asymmetryroot 'data/' args.relevance_dir '/' data_type '/'];
end

%Specify directory to save organisation maps to
radial_dir = [asymmetryroot 'data/' args.radial_dir '/' data_type '/'];
if ~exist(radial_dir, 'dir')
    mkdir(radial_dir);
end

%Work out number of regions and number to inlcude in each job
num_regions = length(line_list);
regions_per_job = ceil(num_regions / args.num_jobs);

%Work out start and end indices for this job
start_idx = (args.task_id - 1)*regions_per_job + 1;
end_idx = min(num_regions, args.task_id*regions_per_job);

display(['Computing maps for images ' num2str(start_idx) ' to ' num2str(end_idx)]);
for ii = start_idx:end_idx;
    
    %Try and load line and orientation maps
    try 
        line_prob = load_uint8([line_dir line_list(ii).name]);
    catch
        display(['problem loading line map ' line_list(ii).name]);
        continue;
    end
    ori_name = dir([ori_dir '*' mam_names{ii} '*.mat']);
    try 
        line_ori = load_uint8([ori_dir ori_name(1).name]);
    catch
        display(['problem loading orientation map ' line_list(ii).name]);
        continue;
    end
    
    if resize
        %Resize the data maps
        line_prob = imresize(line_prob, resize, 'bilinear');
        line_ori = imresize(line_ori, resize, 'bilinear');
    end
    
    %if a relevance maps dir is specified, load in the map and sclae the
    %line probabilities accordingly
    if ~isempty(args.relevance_dir)
        relevance_name = dir([relevance_dir '*' mam_names{ii} '*.mat']);
        try 
            line_relevance = load_uint8([relevance_dir relevance_name(1).name]);
        catch
            display(['problem loading relevance map ' line_list(ii).name]);
            continue;
        end
        line_relevance = imresize(line_relevance, size(line_prob), 'bilinear');
        line_prob = line_prob .* line_relevance;
    end
    
    %Compute organisation map for varying distance weightings
    for jj = 1:length(distance_ranges)
        g_width = distance_ranges(jj);
        
        if g_width == inf
            fname = [radial_dir line_list(ii).name(1:6) '_rad_map_inf.mat'];
            angle_bands = radial_line_projection(line_prob, line_ori, [angular_resolution angular_bands]);

        else
            fname = [radial_dir line_list(ii).name(1:6) '_rad_map_' zerostr(g_width,3) '.mat'];
            display(['Computing organisation map ' fname]);
            [angle_bands] = radial_line_projection(...
                line_prob, line_ori, [angular_resolution angular_bands],...
                fspecial('gaussian', [1 5*g_width], g_width));
        end
        save_uint8(fname, angle_bands);
    end
end
