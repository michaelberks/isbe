function [] = compute_mass_template_map_batch(data_type, varargin)
%COMPUTE_MASS_TEMPLATE_MAP_BATCH *Insert a one line summary here*
%   [] = compute_mass_template_map_batch()
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
% Created: 22-Oct-2010
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester

if ispc
    %task and job identifiers
    task_id = 1;
    num_jobs = 1;
    view = [];
    
    %location of data maps
    mammo_dir = 'mammograms';
    mask_dir = 'masks';
    template_dir = 'template_maps';

    %arguments controlling computation of radial maps
    radii_scales = 5*(6:18);
    radii_ratio = 1.3;
    resize = 0.5;
    
else
%Set defaults for unix system - environment variables listed in .bashrc (and utils/mb_bash.m)
    
    %task and job identifiers
    [z task_id] = unix('echo $SGE_TASK_ID'); task_id = str2num(task_id);
    [z num_jobs] = unix('echo $NUM_JOBS'); num_jobs = str2num(num_jobs);
    [z view] = unix('echo $VIEW'); view(end) = [];
    
    %location of data maps
    [z mammo_dir] = unix('echo $MAMMO_DIR'); mammo_dir(end) = [];
    [z mask_dir] = unix('echo $MASK_DIR'); mask_dir(end) = [];
    [z template_dir] = unix('echo $TEMPLATE_DIR'); template_dir(end) = [];
    
    %arguments controlling computation of radial maps
    [z radii_scales] = unix('echo $RADII_SCALES'); radii_scales = str2num(radii_scales);
    [z radii_ratio] = unix('echo $RADII_RATIO'); radii_ratio = str2num(radii_ratio);
    [z resize] = unix('echo $RESIZE'); resize = str2num(resize);
    
    clear z;
end

args = u_packargs(varargin,... % the user's input
    '0', ... % non-strict mode
    'task_id', task_id, ...
    'num_jobs', num_jobs, ...
    'view', view,...
    'mammo_dir', mammo_dir, ...
    'mask_dir', mask_dir, ...
    'template_dir', template_dir,...
    'radii_scales', radii_scales,...
    'radii_ratio', radii_ratio,...
    'resize', resize);
    
%Get list of line maps
mammo_dir = [asymmetryroot 'data/' args.mammo_dir '/' data_type '/'];
mammo_list = dir([mammo_dir '*' args.view '*.mat']);
mam_names = get_mammo_info(mammo_list);

%Get list of masks (if included)
if ~isempty(args.mask_dir)
    mask_dir = [asymmetryroot 'data/' args.mask_dir '/' data_type '/'];
end

%Specify directory to save organisation maps to
template_dir = [asymmetryroot 'data/' args.template_dir '/' data_type '/'];
if ~exist(template_dir, 'dir')
    mkdir(template_dir);
    mkdir([template_dir '/scales']);
end

%Work out number of regions and number to inlcude in each job
num_regions = length(mammo_list);
regions_per_job = ceil(num_regions / args.num_jobs);

%Work out start and end indices for this job
start_idx = (args.task_id - 1)*regions_per_job + 1;
end_idx = min(num_regions, args.task_id*regions_per_job);

display(['Computing maps for images ' num2str(start_idx) ' to ' num2str(end_idx)]);
for ii = start_idx:end_idx;
    display(['Computing templates for ' mam_names{ii} ' (' num2str(ii) ')']);
    
    %Load mammogram
    mammo = u_load([mammo_dir mammo_list(ii).name]);
    
    %Load mask
    if ~isempty(args.mask_dir)
        mask_name = dir([mask_dir '*' mam_names{ii} '*.mat']); 
        mask = u_load([mask_dir mask_name(1).name]);    
    end
    
    if args.resize
        %Resize the data maps
        mammo = imresize(mammo, args.resize, 'bilinear');
        if ~isempty(args.mask_dir)
            mask = imresize(mask, args.resize, 'bilinear');
        end
    end
    
    %Compute template map and associated scales
    [template_map scales_map] = compute_mass_template_map(mammo, args.radii_scales, args.radii_ratio);
    template_map(~mask) = 0; %#ok
    scales_map(~mask) = 0; %#ok
    
    %Save maps
    save([template_dir mam_names{ii} '_template.mat'], 'template_map');
    save([template_dir '/scales/' mam_names{ii} '_scales.mat'], 'scales_map');

end
display('All mammograms successfully processed');
