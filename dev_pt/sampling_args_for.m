function [sampling_args_out] = sampling_args_for(decomp_type,sampling_args_in)
%SAMPLING_ARGS_FOR Returns subset of sampling_args_in specific to
%decomp_type
%   [sampling_args_out] = sampling_args_for(decomp_type,sampling_args_in)
%
% Inputs:
%      decomp_type - Decomposition type (e.g. 'dt','g2d')
%
%      sampling_args_in - Input structure
%
%
% Outputs:
%      sampling_args_out - Subset of sampling_args_in specific to
%      decomp_type
%
%
% Example:
%      sampling_args = sampling_args_for('dt',sampling_args);
%
% Notes:
%
% See also:
%
% Created: 20-Apr-2011
% Author: Phil Tresadern 
% Email : philip.tresadern@manchester.ac.uk 
% Phone : +44 (0)161 275 5114 
% Copyright: (C) University of Manchester 

switch args.decomp_type
    case 'dt'
		sampling_args_out
        if length(args.num_levels) == 1
            sampling_args.levels = 1:sampling_args_in.num_levels;
        else
            sampling_args.levels = sampling_args_in.num_levels;
        end
        sampling_args.feature_shape = sampling_args_in.feature_shape;
        sampling_args.feature_type = sampling_args_in.feature_type;
        sampling_args.do_max = args.do_max;
        sampling_args.rotate = args.rotate;
        sampling_args.win_size = args.win_size;
    case 'mono'
        sampling_args.win_size = args.win_size;
    case {'g2d', 'g2di', 'clover', 'haar'}
        sampling_args.win_size = args.win_size;
    case 'linop'
        sampling_args.win_size = args.win_size;
        sampling_args.num_levels = args.num_levels;
        sampling_args.num_angles = args.num_angles;
        sampling_args.do_max = args.do_max;
        sampling_args.rotate = args.rotate;
    case 'pixel'
        sampling_args.win_size = args.win_size;
    otherwise
        warning(['decomposition type: ', args.decomp_type, ' not recognised, using DT-CWT']); %#ok
        sampling_args.feature_shape = args.feature_shape;
        sampling_args.feature_type = args.feature_type;
        sampling_args.do_max = args.do_max;
        sampling_args.rotate = args.rotate;
        sampling_args.win_size = args.win_size;
end
