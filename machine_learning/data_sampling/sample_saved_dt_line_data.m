function [X y] = sample_saved_dt_line_data(varargin)
%SAMPLE_SAVED_DT_LINE_DATA *Insert a one line summary here*
%   [] = sample_saved_dt_line_data(varargin)
%
% SAMPLE_SAVED_DT_LINE_DATA uses the U_PACKARGS interface function
% and so arguments to the function may be passed as name-value pairs
% or as a struct with fields with corresponding names. The field names
% are defined as:
%
% Mandatory Arguments:
%       data_dir: full path to the directory containg the saved DT data
%
% Optional Arguments:
%
% Outputs:
%       X: DT samples with representation as specified by the
%       option arguments
%
%       y: Saved outputs (e.g. line class or orientation) associated with each sample
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 21-Feb-2011
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

% Unpack the arguments:
args = u_packargs(varargin, '0',...
    {'saved_data_dir',...
    'num_samples',...
    'task_id'},...
	'pts_per_image', [], ...
    'id_offset', 0,...
    'feature_shape', 'rect',...
    'feature_type', 'conj',...
    'num_levels', 5,...
    'do_max', 0,...
    'rotate', 0,...
    'win_size', 3,...
    'pca', []);
clear varargin;

%Load the training data
args.saved_data_dir = [args.saved_data_dir '/' zerostr(args.task_id+args.id_offset,3) '/'];

if length(args.num_levels) == 1
    num_levels = args.num_levels;
    levels = 1:args.num_levels;
    
else
    num_levels = length(args.num_levels);
    levels = args.num_levels;
end
%Workout dimensions on data based on specified DT representation
if ~isempty(args.pca);
    %pca should have fields mean and modes
    size_sample_vector = size(args.pca.modes,2);
else
    size_sample_vector = num_levels*compute_dt_feature_size(args);
end

if isempty(args.pts_per_image)
	% Workout how many blocks of 10k we need to fill num_samples
	num_blocks = ceil(args.num_samples / 1e4);

	% Pre-allocate data
	X = zeros(num_blocks*1e4, size_sample_vector);
	y = zeros(num_blocks*1e4, 1);
	
	% copy a whole block at a time
	pts_per_block = 1e4;
else
	% Pre-allocate data
	X = zeros(args.num_samples, size_sample_vector);
	y = zeros(args.num_samples, 1);

	% Workout how many blocks of 10k we need to fill num_samples as best we can
	load([args.saved_data_dir,'parameters.mat']);
	num_images = args.num_samples/args.pts_per_image;
	imgs_per_block = length(parameters)/20;
	num_blocks = min(ceil(num_images/imgs_per_block),20);
	
	% randomly select pts_per_block points from each block
	pts_per_block = ceil(args.num_samples/num_blocks);
end

%Loop through blocks, loading data and converting to the specified
%representation
dest_rows = 1:pts_per_block; % rows to copy to
for ii = 1:num_blocks
    
    %load data blocks
    X_10k = u_load([args.saved_data_dir 'X_' zerostr(ii,2) '.mat']);
    y_10k = u_load([args.saved_data_dir 'y_' zerostr(ii,2) '.mat']);
    
	if ~isempty(args.pts_per_image)
		% choose rows from which to pick samples
		src_rows = randperm(size(X_10k,1));
		src_rows = src_rows(1:length(dest_rows));
	
		% reduce to randomly selected rows
		X_10k = X_10k(src_rows,:,:,:);
		y_10k = y_10k(src_rows,:);
	end
	
    %Discard unwanted window dimensions if using 1x1 windows
    if args.win_size == 1
        X_10k = X_10k(:,5,:,:);
    end

    %Discard unwanted levels
    X_10k = X_10k(:,:,:,levels);

    %Convert DT representation and allocate to main data
    rep_args = get_substructure(args,...
        {'feature_type', 'do_max', 'win_size', 'pca'});
    X(dest_rows,:) = convert_dt_representation(X_10k, rep_args);
    
    %Copy y values
    y(dest_rows,:) = y_10k;
	
	% update destination rows
	dest_rows = dest_rows + pts_per_block;
	dest_rows(dest_rows>args.num_samples) = [];
	
	% if we're done then break
	if isempty(dest_rows), break; end
end

%Throwaway additional data if a small number has been specified
if size(X,1)>args.num_samples
    X(args.num_samples+1:end,:) = [];
    y(args.num_samples+1:end) = [];
end