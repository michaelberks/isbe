function [] = generate_bar_training_data(varargin)
%GENERATE_BAR_TRAINING_DATA *Insert a one line summary here*
%   [] = generate_bar_training_data(varargin)
%
% GENERATE_BAR_TRAINING_DATA uses the U_PACKARGS interface function
% and so arguments to the function may be passed as name-value pairs
% or as a struct with fields with corresponding names. The field names
% are defined as:
%
% Mandatory Arguments: Chen test comment
%
% Optional Arguments:
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 29-Jan-2010
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester 

% Unpack the arguments:
args = u_packargs(varargin,... % the user's input
			 '0', ... % non-strict mode
			 {'num_images',... % the mandatory arguments
             'image_dir'}, ...
             'num_levels', 4,...
             'compute_dt', 0);
         
if ~isempty(args.image_dir) && ~strcmp(args.image_dir(end), filesep)
    args.image_dir = [args.image_dir filesep];
end
mkdir(args.image_dir);

for ii = 1:args.num_images
    
    %Generate bar image and label
    [image_out, label, parameters] = create_bar_image; %#ok
    
    max_dt = []; %#ok
    if args.compute_dt
        % Create DT-CWT of image
        dt = dtwavexfm2(image_out, args.num_levels, 'near_sym_b', 'qshift_b');

        %interpolation DT-CWT to full pixel grid
        dt_inter = dt_to_full_image(dt); clear dt;

        %get the maximum response across orientations
        max_dt = squeeze(max(dt_inter, [], 3)); clear dt; %#ok
    end
    
    %save the image data
    save([args.image_dir 'bar', zerostr(ii,3)], 'image_out', 'label', 'parameters', 'max_dt');
    display(['Saved image ', num2str(ii)]);
end