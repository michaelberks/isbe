function [num_channels, channels] = get_channels(args)
% Given sampling arguments, return the number of colour channels used and
% their indices

switch args.rgb_channel
	case 'all'
		%num_channels = 3;
		%channels = 1:3;     
		num_channels = 2;
		channels = 1:2;
		
	case 'r'
		num_channels = 1;
		channels = 1;
		
	case 'g'
		num_channels = 1;
		channels = 2;
		
	case 'b'
		num_channels = 1;
		channels = 3;
		
	case 'rgb'
        % image flattened using rgb2gray()
		num_channels = 1;
        channels = 0;
		
	otherwise
		error(['RGB channel ' args.rgb_channel ' not recognised.']);
end
