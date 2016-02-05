function [args_out] = get_output_args_from(args_in, args_out)
% Given a set of arguments, return only those arguments that are 
% relevant to the predicted outputs (or append them if args_out is supplied)

% This could be rewritten to return a structure with a labels_dir field or
% similar rather than distinct orientation/width fields.

% Deal with undefined inputs
if ~exist('args_out','var'), args_out = []; end

% Arguments common to all output types
args_out.output_type = args_in.output_type;

switch args_in.output_type
    case {'detection', 'centre_detection'}  
        args_out.bg_ratio = args_in.bg_ratio;

    case {'orientation', 'centre_orientation', ...
          'mixed_orientation', 'mixed_centre_orientation'}
        % 100% positive samples 
        args_out.bg_ratio = args_in.bg_ratio;
        args_out.ori_dir = ...
            prettypath([args_in.image_root '/' args_in.ori_dir]);
        args_out.width_dir = [];

    case {'width'}
        % 100% positive samples 
        args_out.bg_ratio = 0;
        args_out.ori_dir = [];
        args_out.width_dir = ...
            prettypath([args_in.image_root '/' args_in.width_dir]);
        
    case {'class_label'}
        args_out.bg_ratio = args_in.bg_ratio;
        args_out.ori_dir = [];
        args_out.width_dir = [];
        args_out.class_label_dir = ...
            prettypath([args_in.image_root '/' args_in.class_label_dir]);
        
    otherwise
        error(['Unknown prediction type:', args_in.output_type]);
end

