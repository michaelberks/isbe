function samples_per_channel = get_samples_per_channel(decomposition_args)
% Given a set of sampling arguments, return the number of samples per
% colour channel

% This might actually be better integrated with get_decomposition_args_from()

% NB: Note the use of length(args.levels) and args.num_levels - they don't
% always correspond for historical reasons :(

% NB: This does not use the do_max argument to reduce the n_levels to 1

% If PCA is used then just return the number of modes and we're done
if isfield(decomposition_args, 'pca') && ~isempty(decomposition_args.pca)
    samples_per_channel = size(decomposition_args.pca,2);
    return;
end

% check decomp_type is valid - note this will also convert any old format
% strings to the new cellstr format
decomposition_args.decomp_type = ...
    check_decomp_type(decomposition_args.decomp_type);

%We can now choose any set of decompositions to combine by having
%decomp_type as a cell array of strings - in which case the samples per
%channel is equal to the sum of each decomp type
if length(decomposition_args.decomp_type) > 1
    samples_per_channel = 0;
    for ii = 1:length(decomposition_args.decomp_type)
        d_args = decomposition_args;
        d_args.decomp_type = d_args.decomp_type(ii);
        samples_per_channel = samples_per_channel + ...
            get_samples_per_channel(d_args);
    end
    return;
end

%Otherwise decomp_type is a string for which we need to compute the samples per channel

%Now uses separate function to compute the number samples in each pixel,
%level and orientation subband block
[n_bands n_pixels n_levels] = get_sample_block_sizes(decomposition_args);

%Samples per channel is just the product of these
samples_per_channel = n_bands * n_pixels * n_levels;

%--------------------------------------------------------------------------
function samples_per_channel = old_func(decomposition_args) %#ok

if ~isempty(decomposition_args.pca)
    samples_per_channel = size(decomposition_args.pca,2);
    return;
end

% check decomp_type is valid
check_decomp_type(decomposition_args.decomp_type);

% number of pixels in window
n_pixels = decomposition_args.win_size^2;

% Get dual tree feature size
switch decomposition_args.decomp_type
	case {'dt', 'dtg2'}
        n_levels = length(decomposition_args.levels);
		dt_feature_size = compute_dt_feature_size(decomposition_args);
end

% Get number of sigmas
switch decomposition_args.decomp_type
	case {'dtg2', 'g1d', 'g1da', 'g12d', 'g2dh', 'g2d', 'clover', 'haar', 'g2da', 'g2dg'...
          'gabor', 'gaborg'}
        n_sigmas = length(decomposition_args.sigma_range);
end

% arguments specific to decomp_types
switch decomposition_args.decomp_type
	case 'dt'
		samples_per_channel = n_levels * dt_feature_size;

	case 'dtg2'
		samples_per_channel = n_levels * dt_feature_size + ...
				  			  3 * n_sigmas * n_pixels;
		
	case 'mono'
        n_levels = decomposition_args.num_levels;
        samples_per_channel = 3 * n_levels * n_pixels;

	case {'g2d', 'clover'}
		samples_per_channel = 3 * n_sigmas * n_pixels;

	case {'haar', 'g1d'}
		samples_per_channel = 2 * n_sigmas * n_pixels;

	case 'g12d'
		samples_per_channel = 5 * n_sigmas * n_pixels;

	case 'g2dh'
		samples_per_channel = 7 * n_sigmas * n_pixels;
        
    case 'g2dg'
		samples_per_channel = 4 * n_sigmas * n_pixels;
        
    case 'g2da'
		samples_per_channel = decomposition_args.num_angles * n_sigmas * n_pixels;
        
	case 'g2di'
		samples_per_channel = 3 * decomposition_args.sigma_range(2) * n_pixels;

	case 'linop'
        if decomposition_args.do_max
			n_angles = 1;
        else
			n_angles = decomposition_args.num_angles;
        end
        n_levels = decomposition_args.num_levels;
		samples_per_channel = n_angles * n_levels * n_pixels;

    case 'gabor'
        % Real and imaginary parts
		samples_per_channel = decomposition_args.num_angles * n_sigmas * ...
                              n_pixels * 2;
                          
    case 'gaborg'
        % Real and imaginary parts
		samples_per_channel = (decomposition_args.num_angles * n_sigmas * n_pixels * 2) + ...
                              (n_sigmas * n_pixels);

    case 'pixel'
		samples_per_channel = n_pixels;
		
	otherwise,
        error(['samples_per_channel not defined for ',decomposition_args.decomp_type]);
end