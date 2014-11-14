function [args_out] = get_decomposition_args_from(args_in, args_out)
% Given a set of arguments, return only those arguments that are 
% relevant to the decomposition type (or append them if args_out is supplied)

% deal with undefined inputs
if ~exist('args_out','var'), args_out = []; end
       
% check decomp_type is valid - this also converts old string types to the
% new cellstr form
args_out.decomp_type = check_decomp_type(args_in.decomp_type);

% This is more like a sampling parameter and should be written out over
% time - MB disagree
args_out.win_size = args_in.win_size;

args_out.rgb_channel = args_in.rgb_channel;
args_out.normalise = args_in.normalise;

% arguments specific to decomp_types - note we're not going to check that
% specific combinations of decomp types don't overwrite each others args in a bad way,
% we'll just have to trust the user has asked for something sensible
for ii = 1:length(args_out.decomp_type)
    switch args_out.decomp_type{ii}
        case 'dt'
            args_out.num_levels = args_in.num_levels;
            if length(args_in.num_levels) == 1
                args_out.levels = 1:args_in.num_levels;
            else
                args_out.levels = args_in.num_levels;
            end
            args_out.feature_shape = args_in.feature_shape;
            args_out.feature_type = args_in.feature_type;
            args_out.do_max = args_in.do_max;
            args_out.rotate = args_in.rotate;
            args_out.use_nag = args_in.use_nag;

        case 'linop',
            args_out.num_levels = args_in.num_levels;
            args_out.num_angles = args_in.num_angles;
            args_out.do_max = args_in.do_max;
            args_out.rotate = args_in.rotate;       

        case {'gabor', 'gabori'}
            args_out.num_angles = args_in.num_angles;
            args_out.sigma_range = args_in.sigma_range;	
            args_out.do_max = args_in.do_max;
            args_out.rotate = args_in.rotate;
            args_out.feature_type = args_in.feature_type;

        case 'mono',
            args_out.num_levels = args_in.num_levels;
            args_out.min_wavelength = args_in.min_wavelength;
            args_out.onf = args_in.onf;

        case {'g1d', 'g2d', 'g2di', 'h2d', 'g', 'haar', 'h2di'},
            args_out.sigma_range = args_in.sigma_range;	

        case {'g2da', 'h2da', 'g2dia', 'h2dia'}
            args_out.sigma_range = args_in.sigma_range;
            args_out.num_angles = args_in.num_angles;
            args_out.do_max = args_in.do_max;
            args_out.rotate = args_in.rotate;

        case {'pixel'}
            args_out.subtract_mean = args_in.subtract_mean;
    end
end

% load in saved PCA data if necessary
if ~isempty(args_in.pca_filename)
    % Not sure about using image_root here but it hardly gets used anyway
    pca_path = prettypath([args_in.image_root '/' args.pca_filename]);
    args_out.pca = u_load(pca_path);
else
    args_out.pca = [];
end
