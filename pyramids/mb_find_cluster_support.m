function [total_hist] = mb_find_cluster_support(varargin)
%
% MB_FIND_CLUSTER_SUPPORT - go through set of pyramids and find counts for
% each cluster from GMM model in given band of pyramid
%
% This function performs the steerable pyramid decomposition
% on a set of shape normalised breasts, and saves the pyramids in a cell 
% array form that is easier to work with than the Simoncelli code's
% vector-based representation.
%
% Mandatory Arguments:
%
%	'PyramidDir'
%		- the directory containing the pyramids.
%
%   'Level'
%	'Orientation'
%   'PyramidModel'
%   'Border'
%
% Return Values:
%   'cluster_hist'
%
% None.

% unpack the arguments
args = u_packargs(varargin, 0, ...
			{'PyramidDir',...
            'Level',...
            'Orientation',...
            'PyramidModel',...
            'Border'},...
            'Plot', 0,...
            'SaveDir', []);
clear varargin

%check pyrmaid directory is filesep terminated
if ~strcmp(args.PyramidDir(end), filesep)
	args.PyramidDir = [args.PyramidDir filesep];
end

%get number of clusters in model
num_clust = size(args.PyramidModel.Means, 1);

% get a listing of the input files
input_list = dir([args.PyramidDir '*pyramid*']);

%pre-allocate hist counts
total_hist = zeros(1, num_clust);

% loop over each pyramid, select the given band and assign nearest cluster
for ii = 1 : length(input_list)
    %load pyramid
    pyramid = u_load([args.PyramidDir, input_list(ii).name]);
    
    %extract band and clear rest of pyramid from memory
    pyr = pyramid{args.Level, args.Orientation};
    clear pyramid
    
    %Get rid of the edges of the image we're not interested in
    pyr = pyr(args.Border+1:end-args.Border, args.Border+1:end-args.Border);
    
    % Assign cluster from model to pyr
    [cluster_image] = assign_cluster_to_image(args.PyramidModel, pyr);
    
    %histogram count the cluster in pyr
    curr_hist = hist(cluster_image(:), 0:num_clust); %0 accounts for edges we can't calculate
    total_hist = total_hist + curr_hist(2:end);
    
    if ~isempty(args.SaveDir)
        c_name = [args.SaveDir, input_list(ii).name, '_', args.Level, '_', args.Orientation];
        save(c_name, 'cluster_image'); 
    end
    
    display(['Completed assigning to image ', num2str(ii)]);
    
    if args.Plot
        figure; imagesc(cluster_image); axis image;
    end
end
    

    
    
    
    
    
    
