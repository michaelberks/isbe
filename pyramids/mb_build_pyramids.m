function [] = mb_build_pyramids(varargin)
%
% MB_BUILD_PYRAMIDS - perform steerable pyramid decomposition
% on a set of images.
%
% This function performs the steerable pyramid decomposition
% on a set of shape normalised breasts, and saves the pyramids in a cell 
% array form that is easier to work with than the Simoncelli code's
% vector-based representation.
%
% Mandatory Arguments:
%
%	'ImageDir'
%		- the directory containing the images.
%
% 'OutputDir'
%		- the directory to save the steerable pyramids in.
%
% Optional Arguments:
%
%	'NumPyramidLevels'
%		- the number of pyramid levels. Defaults to 5.
%
% 'NumOrientations'
%		- the number of orientations. Defaults to 5.
%
% 'ImageFormat'
%		- the format of the input images, defaults to BMP
%
% Return Values:
%
% None.

% unpack the arguments
args = u_packargs(varargin, 0, ...
			{'ImageDir', 'OutputDir'}, ... % Mandatory arguments
			'NumLevels', 5, ...
			'NumOrientations', 5, ...
			'ImageFormat', 'bmp' ...
			);

%Make the output directory if it doesn't already exist
if ~isdir(args.OutputDir)
	mkdir(args.OutputDir);
end

% ensure the directory names are filesep-terminated
if ~strcmp(args.OutputDir(end), filesep)
	args.OutputDir = [args.OutputDir filesep];
end

if ~strcmp(args.ImageDir(end), filesep)
	args.ImageDir = [args.ImageDir filesep];
end

% get a listing of the input files
input_list = dir([args.ImageDir '*.' args.ImageFormat]);

% loop over each image, performing the decomposition
for i = 1 : length(input_list)
	% load the image
	image = double(imread([args.ImageDir input_list(i).name]));
	
	% perform the decomposition
	[pyr, pind] = mb_buildSFpyr(image, args.NumLevels, args.NumOrientations-1);
	
	% clear image, as it's big
	clear('image');
	
    % now convert the pyramids to the CJR format
    pyr = mb_change_pyramid_form(pyr, pind); %#ok
    
	% save the decomposition
	save([args.OutputDir input_list(i).name(1:end-length(args.ImageFormat)-1),...
        '_pyramid.mat'], 'pyr');
	
	% clear the pyramid, as it's big
	clear('pyr');
	clear('pinds')
	
end
