function [] = mb_build_dual_tree(varargin)
%
% MB_BUILD_DUAL_TREE - perform steerable pyramid decomposition
% on a set of images.
%
% This function performs the dual-tree discrete wavelet decomposition
% on a set of images
%
% Mandatory Arguments:
%
%	'ImageDir'
%		- the directory containing the images.
%
% 'OutputDir'
%		- the directory to save the decomposition in.
%
% Optional Arguments:
%
%	'NumLevels'
%		- the number of pyramid levels. Defaults to 5.
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
			'NumLevels', 6, ...
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

%decide if we loading a mat file or image file
mat_file = false;
if strcmp(args.ImageFormat, 'mat')
    mat_file = true;
end

% get a listing of the input files
input_list = dir([args.ImageDir '*.' args.ImageFormat]);

% loop over each image, performing the decomposition
for i = 1 : length(input_list)
	
    % load the image
    if mat_file
        image_in = double(u_load([args.ImageDir input_list(i).name]));
    else
        image_in = double(imread([args.ImageDir input_list(i).name]));
    end
	
	% perform the decomposition
	dual_tree = dtwavexfm2(image_in, args.NumLevels, 'near_sym_b','qshift_b'); %#ok
	
	% clear image, as it's big
	clear('image_in');
    
	% save the decomposition
	save([args.OutputDir input_list(i).name(1:end-length(args.ImageFormat)-1),...
        '_dual_tree.mat'], 'dual_tree');
	
	% clear the decomposition, as it's big
	clear('dual_tree');
	
end
