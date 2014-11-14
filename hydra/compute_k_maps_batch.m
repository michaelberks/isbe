function [] = compute_k_maps_batch(data_type, method_type, varargin)
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
warning('off', 'load_uint8:missing_variables');

args = u_packargs(varargin,... % the user's input
    '0', ... % non-strict mode
    'task_id',				unixenv('SGE_TASK_ID',1), ...
    'num_jobs',             unixenv('NUM_JOBS',1), ...
    'view',                 unixenv('VIEW',[]), ...
    'mam_names',            unixenv('MAM_NAMES',[]), ...
    'line_dir',             unixenv('LINE_DIR',[]), ...
    'ori_dir',              unixenv('ORI_DIR', 'orientation_maps'), ...
    'mask_dir',             unixenv('MASK_DIR','masks'), ...
    'pectoral_dir',         unixenv('PECTORAL_DIR','masks_pectoral'), ...
    'relevance_dir',        unixenv('RELEVANCE_DIR',[]), ...
    'radial_dir',           unixenv('RADIAL_DIR','k_stellate_maps'), ...
    'num_angles',           unixenv('ANGULAR_RESOLUTION', 24), ...
    'spacing',              unixenv('SPACING', 2), ...
    'offset_x',             unixenv('OFFSET_X', 0), ...
    'offset_y',             unixenv('OFFSET_Y', 0), ...
    'thresh',               unixenv('THRESH', []), ...
    'px_per_mm',            unixenv('PX_PER_MM', 100/9),...
    'sigma_range',          unixenv('SIGMA_RANGE', [3.2 4.0 5.4 6.2 7.8]),...
    'r_min',                unixenv('R_MIN', 4),...
    'R',                    unixenv('R_KARSSEMEIJER', 2),...
    'all_pixels',           unixenv('ALL_PIXELS_K', 0),...
    'resize',               unixenv('RESIZE',[]));

display(args);

%Get k-maps parameters from arguments
num_angles = args.num_angles;
spacing = args.spacing;
offset_x = args.offset_x;
offset_y = args.offset_y;
r_max = 3*round(args.px_per_mm*args.sigma_range);
r_min = args.px_per_mm*args.r_min;
R = args.px_per_mm*args.R;

%Get orientation dir from input arguments
ori_dir = [asymmetryroot 'data/' args.ori_dir '/' method_type '/' data_type '/'];

%Either use prescribed mammogram names, or get listing of the mammograms in
%this dir
if isempty(args.mam_names)
    mam_names = get_mammo_info(dir([ori_dir '*' args.view '*.mat']));
else
    mam_names = u_load([asymmetryroot 'data/' args.mam_names]);
end

%Get names of orientation maps
ori_names = match_mammo_names(ori_dir, mam_names);


%check if we're using separate line maps and if so, get their names
if ~isempty(args.line_dir)
    line_dir = [asymmetryroot 'data/' args.line_dir '/' method_type '/' data_type '/'];
    line_names = match_mammo_names(line_dir, mam_names);
else
    line_names = [];
end
    
%Check we're using masks and if so, get their names
if ~isempty(args.mask_dir)
    mask_dir = [asymmetryroot 'data/' args.mask_dir '/' data_type '/'];
    mask_names = match_mammo_names(mask_dir, mam_names);
    
    %Check we're using pectoral masks and if so, get their names
    if ~isempty(args.pectoral_dir)
        pectoral_dir = [asymmetryroot 'data/' args.pectoral_dir '/' data_type '/'];
        pectoral_names = match_mammo_names(pectoral_dir, mam_names);
    else
        pectoral_names = [];
    end
    
else
    mask_names = [];
end

%Check we're using relevance maps and if so, get their names
if ~isempty(args.relevance_dir)
    relevance_dir = [asymmetryroot 'data/' args.relevance_dir '/' data_type '/'];
    relevance_names = match_mammo_names(relevance_dir, mam_names);
else
    relevance_names = [];
end


%Specify directory to save organisation maps to
radial_dir = [asymmetryroot 'data/' args.radial_dir '/' method_type '/' data_type '/'];
if ~exist(radial_dir, 'dir')
    mkdir(radial_dir);
end

%Work out number of regions and number to inlcude in each job
num_regions = length(mam_names);
regions_per_job = ceil(num_regions / args.num_jobs);

%Work out start and end indices for this job
start_idx = (args.task_id - 1)*regions_per_job + 1;
end_idx = min(num_regions, args.task_id*regions_per_job);

display(['Computing maps for images ' num2str(start_idx) ' to ' num2str(end_idx)]);
for ii = start_idx:end_idx;
    
    %Try orientation maps
    try 
        ori_map = load_uint8([ori_dir ori_names{ii}]);
    catch
        display(['problem loading orientation map: "' [ori_dir ori_names{ii}] '"'...
            ' for mammogram ' mam_names{ii}]);
        continue;
    end
    
    %If we're not using separate line maps, compute line map as the
    %magnitude of the orientation map and then take the complex arg of the
    %orientation map
    if isempty(line_names)
        if args.all_pixels
            line_map = true(size(ori_map));
        else
            line_map = abs(ori_map);
        end
    else
        %try loading the line map
        try 
            line_map = load_uint8([line_dir line_names{ii}]);
        catch
            display(['problem loading line map: "' [line_dir line_names{ii}] '"'...
                ' for mammogram ' mam_names{ii}]);
            continue;
        end
        %Check whether line map is complex
        if ~isreal(line_map)
            line_map = abs(line_map);
        end
    end
    
    %Check whether orientation map is complex
    if ~isreal(ori_map)
        ori_map = angle(ori_map);
    end
    
    %Check if we're thresholding the line map to obtain a binary image
    if ~isempty(args.thresh)
        line_map = line_map > args.thresh;
    end
    
    %Check if we've asked to resize the maps
    if args.resize
        %Resize the data maps
        line_map = imresize(line_map, args.resize, 'bilinear');
        ori_map = imresize(ori_map, resize, 'bilinear');
    end
    
    %if a relevance maps dir is specified, load in the map and sclae the
    %line probabilities accordingly
    if ~isempty(relevance_names)
        try 
            relevance_map = load_uint8([relevance_dir relevance_names{ii}]);
        catch
            display(['problem loading relevance map: "' [relevance_dir relevance_names{ii}] '"'...
                ' for mammogram ' mam_names{ii}]);
            continue;
        end
        relevance_map = imresize(relevance_map, size(line_map), 'bilinear');
        line_map = line_map .* relevance_map;
    end
    
    %if mask names are set, load in the mask
    if ~isempty(mask_names)
        try 
            mask = load_uint8([mask_dir mask_names{ii}]);
        catch
            display(['problem loading mask: "' [mask_dir mask_names{ii}] '"'...
                ' for mammogram ' mam_names{ii}]);
            continue;
        end
        mask = imresize(mask, size(line_map));
        
        %Also see if we have a pectoral mask to load in
        if ~isempty(pectoral_names) && ~isempty(pectoral_names{ii})
            try 
                pectoral_mask = load_uint8([pectoral_dir pectoral_names{ii}]);
            catch
                display(['problem loading pectoral mask: "' [pectoral_dir pectoral_names{ii}] '"'...
                    ' for mammogram ' mam_names{ii}]);
                continue;
            end
            pectoral_mask = imresize(pectoral_mask, size(line_map));
            mask = mask & ~pectoral_mask; clear pectoral_mask;
        end
        
        %Discard lines from outside the mask
        line_map(~mask) = 0;
            
    else
        mask = [];
    end
    
    %Compute map of pixels to discard (i.e. near the edge of the image
    %where artefacts cause strong aligned orientations
    [discard_mask] = discard_orientations(ori_map, 'mask', mask);  
    line_map(discard_mask) = 0;  
    
    %Compute the Karssemeijer maps
    [f_i1 f_i2 mask] = karssemeijer_radial_projection_multiscale(...
        line_map, ori_map,...
        'r_min', r_min,...
        'r_max', r_max,...
        'R', R,...
        'num_angles', num_angles,...
        'min_thresh', 20,...
        'spacing', spacing,...
        'offset_x', offset_x,...
        'offset_y', offset_y,...
        'mask', mask);
    
    %Save the maps
    save_uint8([radial_dir mam_names{ii} '_f1.mat'], f_i1);
    save_uint8([radial_dir mam_names{ii} '_f2.mat'], f_i2);
    save([radial_dir mam_names{ii} '_mask.mat'], 'mask');

end
