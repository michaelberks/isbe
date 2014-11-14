function [num_channels, channels] = get_rgb_channels(args_in)
% Given a set of arguments, return the number of colour channels used and
% their indices.

if ~isfield(args,'rgb_channel')
    error('args_in must have rgb_channel field');
end

switch args_in.rgb_channel
    case {'all'}
        num_channels = 3;
        channels = 1:3;
        
    case {'rgb'}
        num_channels = 1;
        % channels not defined here - rgb2gray() is used instead.

    case {'r','mono'}
        num_channels = 1;
        channels = 1;
        
    case {'g'}
        num_channels = 1;
        channels = 2;
        
    case {'b'}
        num_channels = 1;
        channels = 3;
        
    otherwise
        error(['RGB channel identifier ' args_in.rgb_channel ' not recognised.']);
end
