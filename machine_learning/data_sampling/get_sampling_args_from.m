function [args_out] = get_sampling_args_from(args_in, args_out)
% Given a set of arguments, return only those arguments that are 
% relevant to the sampling strategy (or append them if args_out is supplied)

% deal with undefined inputs
if ~exist('args_out','var'), args_out = []; end

% Arguments common to all samplers
args_out.sampling_method = args_in.sampling_method;

args_out.image_root = prettypath(args_in.image_root);
args_out.image_type = args_in.image_type;
args_out.num_samples = args_in.num_samples;
args_out.max_n_images = args_in.max_n_images;
args_out.shift_images = args_in.shift_images;
args_out.shrink_fov = args_in.shrink_fov;
args_out.win_size = args_in.win_size;
args_out.replace_sample = args_in.replace_sample;

% args_out.training_data_dir = ...
%     prettypath([args_in.image_root '/' args_in.training_data_dir]);


switch args_in.sampling_method
    % Generic samplers (ideally the only two we'd need)
    case {'generate_training_data'}
        args_out = get_generate_args(args_in, args_out);
        
    case {'sample_saved_training_data'}
        args_out = get_resample_args(args_in, args_out);
       
    % Mammogram samplers
    case {'sample_vessel_dt_data'}
        bg_root = [asymmetryroot, 'data/synthetic_backgrounds/'];
        args_out.saved_data_dir = [bg_root args_in.bg_dir];
        args_out.id_offset = (args_in.task_id-1) * args_in.n_trees; 
        
		if strcmp(get_username,'ptresadern')
			args_out.pts_per_image = args_in.pts_per_image;
        end
    
    case {'sample_saved_dt_line_data'}
        bg_root = [asymmetryroot, 'data/synthetic_backgrounds/'];
        args_out.saved_data_dir = [bg_root args_in.bg_dir];
        args_out.id_offset = (args_in.task_id-1) * args_in.n_trees; 
        args_out.save_training_data = args_in.save_training_data;

        if strcmp(get_username,'ptresadern')
			args_out.pts_per_image = args_in.pts_per_image;
        end
        
    otherwise,
        error(['Unknown sampling method: ', args_in.sampling_method]);
end


function args_out = get_generate_args(args_in, args_out)
% Return arguments for functions that generate new training data

% deal with undefined inputs
if ~exist('args_out','var'), args_out = []; end

switch(args_in.image_type)
    case {'real'}
        % Sample from real images that are loaded from folders
        args_out.image_dir = ...
            prettypath([args_in.image_root '/' args_in.image_dir]);
        if ~isempty(args_in.fov_mask_dir)
            args_out.fov_mask_dir = ...
                prettypath([args_in.image_root '/' args_in.fov_mask_dir]);
        else
            args_out.fov_mask_dir = [];
        end
        if ~isempty(args_in.fg_mask_dir)
            args_out.fg_mask_dir = ...
                prettypath([args_in.image_root '/' args_in.fg_mask_dir]);
        else
            args_out.fg_mask_dir = [];
        end
        if ~isempty(args_in.probability_dir)
            if args_in.make_resampling_maps
                %Will compute resampling probabilities on the fly from
                %predictions
                args_out.probability_dir = ...
                    prettypath([args_in.image_root ...
                    '/' args_in.prediction_dir '/' args_in.output_type ...
                    '/' args_in.prediction_type '/' args_in.probability_dir]);
            else
                %Will load in existing sampling maps
                args_out.probability_dir = ...
                    prettypath([args_in.image_root '/' args_in.probability_dir]);
            end
            args_out.make_resampling_maps = args_in.make_resampling_maps;
        else
            args_out.probability_dir = [];
        end
        if isfield(args_in, 'precompute_indices')
            args_out.precompute_indices = args_in.precompute_indices;
        else
            args_out.precompute_indices = [];
        end
        
    case {'line', 'grain'}
        % Create synthetic images on the fly
        fields_to_copy = {...
            'pts_per_image'; 'num_bgs';
            'bg_type'; 'bg_size';
            'bg_dir'; 'bg_stem'; 'bg_zeros'; 'bg_fmt'; 'bg_mask_dir'; 
            'orientation_range'; 'width_range'; 'contrast_range'; 
            'decay_rate'; 'line_type'; 'noise_type'; 'noise_params';
        };
        args_out = get_substructure(args_in, fields_to_copy, args_out);
        
    otherwise
        error(['Image type ', args_in.image_type, ' not recognized']);
end


function args_out = get_resample_args(args_in, args_out)
% Return arguments for functions that resample previously generated training 
% data
args_out.training_data = args_in.training_data;
args_out.training_labels = args_in.training_labels;
