function [] = classify_mammo_fold(varargin)
%
%
%
% Inputs:
%
%
% Outputs:
%
% Example: for use on hydra
%   NUM_JOBS=20 FOREST_JOB="'191658'" TEST_IMAGE_DIR="'lines512'" qsub -N class_im -t 1-20 -V matlab_code/trunk/hydra/classify_image_set.sh
%
% Notes:
%
% See also:
%
% Created: 03-Feb-2010
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester
%Set defaults for windows system
if ispc
    
    %task and job identifiers
    task_id = 1;
    fold_id = 1;
    custom_id = [];
    num_folds = 10;
    num_jobs = 10;
    classify_abnormals = 1;
    classify_contra = 1;
    classify_normals = 1;
    
    %arguments controlling data sampling
    use_nag = 1;
    
    abnormal_data = '2004_screening\abnormals';
    normal_data = '2004_screening\normals';
    image_dir = 'mammograms';
    mask_dir = 'masks';
    mammo_rf_dir = 'mammo_rfs';
    
    image_type = '.mat';
    view = [];
    
    %arguments controlling forest classification
    forest_type = 'isbe';
    max_size = 128;
    num_trees = [];
    use_probs = 0;

else
%Set defaults for unix system - environment variables listed in .bashrc (and utils/mb_bash.m)
    %task and job identifiers
    [z task_id] = unix('echo $SGE_TASK_ID'); task_id = str2num(task_id); %#ok
    [z custom_id] = unix('echo $CUSTOM_ID'); custom_id(end) = [];
    [z fold_id] = unix('echo $FOLD_ID'); fold_id = str2num(fold_id); %#ok
    [z num_folds] = unix('echo $NUM_FOLDS'); num_folds = str2num(num_folds); %#ok
    [z num_jobs] = unix('echo $NUM_JOBS'); num_jobs = str2num(num_jobs); %#ok
    [z classify_normals] = unix('echo $CLASSIFY_NORMALS'); classify_normals = str2num(classify_normals); %#ok
    [z classify_abnormals] = unix('echo $CLASSIFY_ABNORMALS'); classify_abnormals = str2num(classify_abnormals); %#ok
    [z classify_contra] = unix('echo $CLASSIFY_CONTRA'); classify_contra = str2num(classify_contra); %#ok
    
    %arguments controliing data sampling
    [z use_nag] = unix('echo $USE_NAG'); use_nag = str2num(use_nag); %#ok
    [z abnormal_data] = unix('echo $ABNORMAL_DATA'); abnormal_data(end) = [];
    [z normal_data] = unix('echo $NORMAL_DATA'); normal_data(end) = [];
    [z image_dir] = unix('echo $IMAGE_DIR'); image_dir(end) = [];
    [z mask_dir] = unix('echo $MASK_DIR'); mask_dir(end) = [];
    [z mammo_rf_dir] = unix('echo $MAMMO_RF_DIR'); mammo_rf_dir(end) = [];
    [z image_type] = unix('echo $IMAGE_TYPE'); image_type(end) = [];
    [z view] = unix('echo $VIEW'); view(end) = [];
    
    %arguments controlling forest classification
    [z forest_type] = unix('echo $FOREST_TYPE'); forest_type(end) = [];
    [z max_size] = unix('echo $MAX_SIZE'); max_size = str2num(max_size); %#ok
    [z num_trees] = unix('echo $NUM_TREES_C'); num_trees = str2num(num_trees); %#ok
    [z use_probs] = unix('echo $USE_PROBS'); use_probs = str2num(use_probs); %#ok

    clear z;
    warning('off', 'load_uint8:missing_variables');
end

%Now use varargin to merge default argument values with user set values
args = u_packargs(varargin,... % the user's input
    '0', ... % non-strict mode
    'task_id', task_id, ...
    'custom_id', custom_id,...
    'fold_id', fold_id, ...
    'num_folds', num_folds, ...
    'num_jobs', num_jobs, ...
    'abnormal_data', abnormal_data,...
    'normal_data', normal_data,...
    'image_dir', image_dir,...
    'mask_dir', mask_dir,...
    'mammo_rf_dir', mammo_rf_dir,...
    'image_type', image_type,...
    'view', view,...
    'forest_type', forest_type, ...
    'max_size', max_size, ...
    'num_trees', num_trees,...
    'use_probs', use_probs,...
    'use_nag', use_nag,...
    'classify_normals', classify_normals,...
    'classify_abnormals', classify_abnormals,...
    'classify_contra', classify_contra);

abnormal_dir = [asymmetryroot, 'data/' args.image_dir '/' args.abnormal_data '/'];
normal_dir = [asymmetryroot, 'data/' args.image_dir '/' args.normal_data '/'];
abnormal_mask_dir = [asymmetryroot, 'data/' args.mask_dir '/' args.abnormal_data '/'];
normal_mask_dir = [asymmetryroot, 'data/' args.mask_dir '/' args.normal_data '/'];

%First find all the abnormal mammograms that match the view type
abnormal_list = dir([abnormal_dir '*' args.view '*' args.image_type]);
abnormal_names = get_mammo_info(abnormal_list);

%Get mathcing list of abnormal masks
abnormal_mask_list = dir([abnormal_mask_dir '*' args.view '*.mat']);

%If there is a meta folder associated with the abnormal dir we must select
%only the mammograms that actually have an abnormality in them, and therefore have meta
%information (as opposed to the contralateral mammogram also stored in the
%abnormal directory)
meta_list = dir([abnormal_dir '/meta/*' args.view '*.mat']);
if ~isempty(meta_list)
    meta_names = get_mammo_info(meta_list);
    [dummy abnormal_idx] = intersect(abnormal_names, meta_names);
    abnormal_list = abnormal_list(abnormal_idx);
    abnormal_mask_list = abnormal_mask_list(abnormal_idx);
end

%Get list of normal mammograms and masks
normal_list = dir([normal_dir '*' args.view '*' args.image_type]);
normal_mask_list = dir([normal_mask_dir '*' args.view '*.mat']);

%Select number of images to use
num_images = min(length(abnormal_list), length(normal_list));

%Workout the number of images to include in ecah fold
fold_size = ceil(num_images / args.num_folds);
start_idx = (args.fold_id-1)*fold_size + 1;
end_idx = min(args.fold_id*fold_size, num_images);

%Now update the mammogram and mask lists
fold_idx = start_idx:end_idx;
abnormal_list = abnormal_list(fold_idx);
normal_list = normal_list(fold_idx);
abnormal_mask_list = abnormal_mask_list(fold_idx);
normal_mask_list = normal_mask_list(fold_idx);

%Workout RF name    
if isempty(args.custom_id)    
    forest_name = [args.view '_' zerostr(args.fold_id,2)];
else
    forest_name = [args.custom_id '_' zerostr(args.fold_id,2)];
end

%Construct paths to images and results directories
abnormal_results_dir = [abnormal_dir '/results/' forest_name(1:end-3) '/'];
normal_results_dir = [normal_dir '/results/' forest_name(1:end-3) '/'];

display(args);
display(abnormal_results_dir);
display(normal_results_dir);
display([asymmetryroot 'data/' args.mammo_rf_dir '/' forest_name '/random_forest.mat']);

%Load in random forest and samplign arguments
forest = u_load([asymmetryroot 'data/' args.mammo_rf_dir '/' forest_name '/random_forest.mat']);
sampling_args = u_load([asymmetryroot 'data/' args.mammo_rf_dir '/' forest_name '/sampling_args.mat']);

if ~exist(abnormal_results_dir, 'dir')
    %Make output directory
    mkdir(abnormal_results_dir);
end
if ~exist(normal_results_dir, 'dir')
    %Make output directory
    mkdir(normal_results_dir);
end

%1. Workout number of images in job
job_size = round(length(abnormal_list) / args.num_jobs);

%2. Workout start and end indices for job
start_i = (args.task_id-1)*job_size + 1;
end_i = min(args.task_id*job_size, length(abnormal_list));

display(['classifying images ' num2str(start_i) ' to ' num2str(end_i)]);

%Increase the recursion limit if it is set below 200
if get(0, 'recursionlimit') < 200
    set(0, 'recursionlimit', 200);
end

%Check for decomposition type
if isfield(sampling_args, 'decomp_type')
    decomp_type = sampling_args.decomp_type;
elseif isfield(sampling_args, 'dist_range')
    decomp_type = 'radial';
else
    decomp_type = 'dt';
end
switch decomp_type
    case 'dt'
        sampling_args_c.num_levels = sampling_args.num_levels;
        sampling_args_c.feature_shape = sampling_args.feature_shape;
        sampling_args_c.feature_type = sampling_args.feature_type;
        sampling_args_c.do_max = sampling_args.do_max;
        sampling_args_c.rotate = sampling_args.rotate;
        sampling_args_c.win_size = sampling_args.win_size;
        sampling_args_c.use_nag = sampling_args.use_nag;
    case 'mono'
        sampling_args_c.num_levels = sampling_args.num_levels;
        sampling_args_c.win_size = sampling_args.win_size;
        sampling_args_c.min_wavelength = sampling_args.min_wavelength;
        sampling_args_c.onf = sampling_args.onf;
    case 'linop'
        sampling_args_c.win_size = sampling_args.win_size;
        sampling_args_c.num_levels = sampling_args.num_levels;
        sampling_args_c.num_angles = sampling_args.num_angles;
        sampling_args_c.do_max = sampling_args.do_max;
        sampling_args_c.rotate = sampling_args.rotate;
    case 'pixel'
        sampling_args_c.win_size = sampling_args.win_size;
    case 'radial'
        sampling_args_c.sigma_range = sampling_args.sigma_range;
        sampling_args_c.angular_res = sampling_args.angular_res;
        num_dists = length(sampling_args.dist_range);
        if isfield(sampling_args, 'do_template')
            sampling_args_c.do_template = sampling_args.do_template;
        else
            sampling_args_c.do_template = 1;
        end
        if isfield(sampling_args, 'do_scale')
            sampling_args_c.do_scale = sampling_args.do_scale;
        else
            sampling_args_c.do_scale = 1;
        end
end
if isfield(sampling_args, 'pca') && ~isempty(sampling_args.pca);
    %pca should have fields mean and modes
    sampling_args_c.pca = u_load(sampling_args.pca);
end

%3. Classify each image using forest
if args.classify_abnormals
    for ii = start_i:end_i

        display(['classifying abnormal mammogram ' num2str(ii)]);

        %load mask
        mask = u_load([abnormal_mask_dir abnormal_mask_list(ii).name]);
        [num breast view] = get_mammo_info(abnormal_list(ii).name);

        if strcmp(decomp_type, 'radial')
            %HACK ALERT - radial maps have been created half the size -
            %this really shouldn't be hard-coded but what the hell...
            mask = imresize(mask, 0.5);
            
            %for radial classification, we only need an input image to get
            %the size in classify_image - so create it as a sparse array
            %the same size as the mask - neat huh?
            mammo = sparse(zeros(size(mask)));
            
            %Need to save map dirs into the sampling args
            sampling_args_c.radial_dir = [sampling_args.radial_dir sampling_args.abnormal_data '/'];
            sampling_args_c.template_dir = [sampling_args.template_dir sampling_args.abnormal_data '/'];
            
            if sampling_args_c.do_template
                %Also need the template map name...
                template_name = dir([sampling_args_c.template_dir '*' num breast view '*.mat']);
                if isempty(template_name)
                    display(['Missing template map for image' num2str(ii) ' - skipping']);
                    continue;
                else
                    sampling_args_c.template_name = template_name(1).name;
                end
            end
            
            if sampling_args_c.do_scale
                %... the scale map name
                scale_name = dir([sampling_args_c.template_dir 'scales/*' num breast view '*.mat']);
                if isempty(scale_name)
                    display(['Missing scale map for image' num2str(ii) ' - skipping']);
                    continue;
                else
                    sampling_args_c.scale_name = scale_name(1).name;
                end
            end
            
            %... and the radial names for each distance range
            radial_names = cell(1, num_dists);
            for jj = 1:num_dists
                rad_name = ...
                    dir([sampling_args_c.radial_dir '*' num breast view '*'...
                        zerostr(sampling_args.dist_range(jj),3), '.mat']);
                if isempty(rad_name)
                    continue;
                else
                    radial_names{jj} = rad_name(1).name;
                end
            end
            if any(cellfun(@isempty, radial_names))
                display(['Missing radial map for image' num2str(ii) ' - skipping']);
                continue;
            else
                sampling_args_c.radial_names = radial_names;
            end
                
            %Don't need to flip right breasts
            flipped = 0;

        else
            %load mammogram
            mammo = u_load([abnormal_dir abnormal_list(ii).name]);
            
            %Check if right breast and flip mask and mammo if necessary
            if strcmpi(breast, 'R')
                mammo = fliplr(mammo);
                mask = fliplr(mask);
            end
            flipped = 1;
        end

        %Classify the image
        [probability_image] = classify_image(...
            'image_in', mammo, ...
            'forest', forest,...
            'sampling_args', sampling_args_c,...
            'forest_type', args.forest_type,...
            'decomp_type', decomp_type,...
            'mask', mask,...
            'num_trees', args.num_trees, ...
            'max_size', args.max_size,...
            'use_probs', args.use_probs);
        
        %Flip back right breasts if necessary
        if flipped && strcmpi(breast, 'R')
            probability_image = fliplr(probability_image);
        end
        
        %Save the image in the output directory
        save_uint8([abnormal_results_dir num breast view  '_class.mat'], probability_image);
    end
end

if args.classify_contra
    for ii = start_i:end_i
        display(['classifying contralateral breast to abnormal mammogram ' num2str(ii)]);

        %load mask
        [num breast view] = get_mammo_info(abnormal_list(ii).name);
        if strcmpi(breast, 'R')
            breast = 'L';
        else
            breast = 'R';
        end
        
        %Also need the template map name...
        mask_name = dir([abnormal_mask_dir '*' num breast view '*.mat']);
        if isempty(mask_name)
            display(['Missing mask map for image' num2str(ii) ' - skipping']);
            continue;
        else
            mask = u_load([abnormal_mask_dir mask_name(1).name]);
        end
            
        if strcmp(decomp_type, 'radial')
            %HACK ALERT - radial maps have been created half the size -
            %this really shouldn't be hard-coded but what the hell...
            mask = imresize(mask, 0.5);
            
            %for radial classification, we only need an input image to get
            %the size in classify_image - so create it as a sparse array
            %the same size as the mask - neat huh?
            mammo = sparse(zeros(size(mask)));
            
            %Need to save map dirs into the sampling args
            sampling_args_c.radial_dir = [sampling_args.radial_dir sampling_args.abnormal_data '/'];
            sampling_args_c.template_dir = [sampling_args.template_dir sampling_args.abnormal_data '/'];
            
            if sampling_args_c.do_template
                %Also need the template map name...
                template_name = dir([sampling_args_c.template_dir '*' num breast view '*.mat']);
                if isempty(template_name)
                    display(['Missing template map for image' num2str(ii) ' - skipping']);
                    continue;
                else
                    sampling_args_c.template_name = template_name(1).name;
                end
            end
            
            if sampling_args_c.do_scale
                %... the scale map name
                scale_name = dir([sampling_args_c.template_dir 'scales/*' num breast view '*.mat']);
                if isempty(scale_name)
                    display(['Missing scale map for image' num2str(ii) ' - skipping']);
                    continue;
                else
                    sampling_args_c.scale_name = scale_name(1).name;
                end
            end
            
            %... and the radial names for each distance range
            radial_names = cell(1, num_dists);
            for jj = 1:num_dists
                rad_name = ...
                    dir([sampling_args_c.radial_dir '*' num breast view '*'...
                        zerostr(sampling_args.dist_range(jj),3), '.mat']);
                if isempty(rad_name)
                    continue;
                else
                    radial_names{jj} = rad_name(1).name;
                end
            end
            if any(cellfun(@isempty, radial_names))
                display(['Missing radial map for image' num2str(ii) ' - skipping']);
                continue;
            else
                sampling_args_c.radial_names = radial_names;
            end
                
            %Don't need to flip right breasts
            flipped = 0;

        else
            %load mammogram
            mammo = u_load([abnormal_dir abnormal_list(ii).name]);
            
            %Check if right breast and flip mask and mammo if necessary
            if strcmpi(breast, 'R')
                mammo = fliplr(mammo);
                mask = fliplr(mask);
            end
            flipped = 1;
        end

        %Classify the image
        [probability_image] = classify_image(...
            'image_in', mammo, ...
            'forest', forest,...
            'sampling_args', sampling_args_c,...
            'forest_type', args.forest_type,...
            'decomp_type', decomp_type,...
            'mask', mask,...
            'num_trees', args.num_trees, ...
            'max_size', args.max_size,...
            'use_probs', args.use_probs);
        
        %Flip back right breasts if necessary
        if flipped && strcmpi(breast, 'R')
            probability_image = fliplr(probability_image);
        end
        
        %Save the image in the output directory
        save_uint8([abnormal_results_dir num breast view  '_class.mat'], probability_image);
    end    
end
%Now repeat the classsification for the normal images
if args.classify_normals
    for ii = start_i:end_i

        display(['classifying normal mammogram ' num2str(ii)]);

        %load mask
        mask = u_load([normal_mask_dir normal_mask_list(ii).name]);
        [num breast view] = get_mammo_info(normal_list(ii).name);

        if strcmp(decomp_type, 'radial')
            %HACK ALERT - radial maps have been created half the size -
            %this really shouldn't be hard-coded but what the hell...
            mask = imresize(mask, 0.5);
            
            %for radial classification, we only need an input image to get
            %the size in classify_image - so create it as a sparse array
            %the same size as the mask - neat huh?
            mammo = sparse(zeros(size(mask)));
            
            %Need to save the map dirs into the sampling args
            sampling_args_c.radial_dir = [sampling_args.radial_dir sampling_args.normal_data '/'];
            sampling_args_c.template_dir = [sampling_args.template_dir sampling_args.normal_data '/'];
            
            if sampling_args_c.do_template
                %Also need the template map name...
                template_name = dir([sampling_args_c.template_dir '*' num breast view '*.mat']);
                if isempty(template_name)
                    display(['Missing template map for image' num2str(ii) ' - skipping']);
                    continue;
                else
                    sampling_args_c.template_name = template_name(1).name;
                end
            end
            
            if sampling_args_c.do_scale
                %... the scale map name
                scale_name = dir([sampling_args_c.template_dir 'scales/*' num breast view '*.mat']);
                if isempty(scale_name)
                    display(['Missing scale map for image' num2str(ii) ' - skipping']);
                    continue;
                else
                    sampling_args_c.scale_name = scale_name(1).name;
                end
            end
            
            %... and the radial names
            radial_names = cell(1, num_dists);
            for jj = 1:num_dists
                rad_name = ...
                    dir([sampling_args_c.radial_dir '*' num breast view '*'...
                        zerostr(sampling_args.dist_range(jj),3), '.mat']);
                if isempty(rad_name)
                    continue;
                else
                    radial_names{jj} = rad_name(1).name;
                end
            end
            if any(cellfun(@isempty, radial_names))
                display(['Missing radial map for image' num2str(ii) ' - skipping']);
                continue;
            else
                sampling_args_c.radial_names = radial_names;
            end
            
            %Don't need to flip right breasts
            flipped = 0;
        else
            %load mammogram
            mammo = u_load([normal_dir normal_list(ii).name]);
            
            %Check if right breast and flip mask and mammo if necessary
            if strcmpi(breast, 'R')
                mammo = fliplr(mammo);
                mask = fliplr(mask);
            end
            flipped = 1;
        end

        %Classify the image
        [probability_image] = classify_image(...
            'image_in', mammo, ...
            'forest', forest,...
            'sampling_args', sampling_args_c,...
            'forest_type', args.forest_type,...
            'decomp_type', decomp_type,...
            'mask', mask,...
            'num_trees', args.num_trees, ...
            'max_size', args.max_size,...
            'use_probs', args.use_probs);
        
        %Flip back right breasts if necessary
        if flipped && strcmpi(breast, 'R')
            probability_image = fliplr(probability_image);
        end
        
        %Save the image in the output directory
        save_uint8([normal_results_dir num breast view  '_class.mat'], probability_image);
    end
end
display('All images successfully classified');
