function [n_bands n_pixels n_levels] = get_sample_block_sizes(decomposition_args)
% Given a set of sampling arguments, return the number of samples per
% colour channel

% This might actually be better integrated with get_decomposition_args_from()

% NB: Note the use of length(args.levels) and args.num_levels - they don't
% always correspond for historical reasons :(

% number of pixels in window
n_pixels = decomposition_args.win_size^2;
switch decomposition_args.decomp_type{1}
    case 'dt'
        n_pixels = compute_dt_feature_size(decomposition_args);
end

% Get number of levels
switch decomposition_args.decomp_type{1}
	case {'g', 'g1d', 'g2d', 'g2da', 'h2d', 'h2da', 'haar', 'gabor'}
        n_levels = length(decomposition_args.sigma_range);
        
    case {'g2di', 'h2di', 'g2dia', 'h2dia', 'gabori'}
        n_levels = decomposition_args.sigma_range(2);
        
    case {'dt'}
        n_levels = length(decomposition_args.levels);
        
    case {'mono','linop'}
        n_levels = decomposition_args.num_levels;
        
    case 'pixel'
        n_levels = 1;
end

%Get number of bands per pixel per level
switch decomposition_args.decomp_type{1}
    case 'dt'
		n_bands = 6;
		
	case {'g2d', 'mono', 'g2di'}
		n_bands = 3 ;

	case {'haar', 'g1d'}
		n_bands = 2;

	case {'h2d', 'h2di'}
		n_bands = 4;
        
    case {'g2da', 'h2da', 'linop', 'gabor', 'gabori', 'g2dia', 'h2dia'}
        n_bands = decomposition_args.num_angles;

    case {'pixel', 'g'}
		n_bands = 1;
end

%Modify n_bands for decomp_types that use do_max
switch decomposition_args.decomp_type{1}
    case {'dt', 'g2da', 'h2da', 'g2dia', 'h2dia', 'linop', 'gabor', 'gabori'}
        if decomposition_args.do_max
            n_bands = 1;
        end
end

%Modify n_bands for complex types
switch decomposition_args.decomp_type{1}
    case {'dt', 'gabor', 'gabori'}
        n_bands = n_bands * compute_complex_feature_size(decomposition_args);
end
