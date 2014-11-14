function [thresh percentiles] = mass_maps_thresh(map_dir, varargin)
%MASS_MAPS_THRESH *Insert a one line summary here*
%   [] = mass_maps_thresh()
%
% Inputs:
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 10-Nov-2010
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester

% Unpack the arguments:
args = u_packargs(varargin,... % the user's input
    '0', ... % non-strict mode
    'mammo_dir', 'mammograms',...
    'data_type', '2004_screening/abnormals',...
    'mammo_list', [],...
    'mask_dir', 'masks',...
    'meta_dir', 'meta',...
    'map_type', [],...
    'dist_ranges', [],...
    'sigma', 0,...
    'percentile', 95,...
    'quiet', 1);
clear varargin;

if args.quiet
    warning('off', 'load_uint8:missing_variables');
end

%form complete paths to map and mammogram directories
map_dir = [asymmetryroot 'data/' map_dir  '/' args.data_type '/'];
mammo_dir = [asymmetryroot 'data/' args.mammo_dir '/' args.data_type '/'];
mask_dir = [asymmetryroot 'data/' args.mask_dir '/' args.data_type '/'];

%get list of mammograms containing abnormalities
if isempty(args.mammo_list)
    mammo_list = dir([mammo_dir '/' args.meta_dir '/*.mat']);
else
    mammo_list = args.mammo_list;
end
mammo_names = get_mammo_info(mammo_list);

num_mammos = length(mammo_names);


if isempty(args.dist_ranges)
    num_dists = 1;
    dist_names{1} = [];
else
    num_dists = length(args.dist_ranges);
    dist_names = cell(num_dists,1);
    for ii = 1:num_dists
        dist_names{ii} = zerostr(args.dist_ranges(ii), 3);
    end
end

percentiles = nan(num_mammos, 2, num_dists);
%
for ii = 1:num_mammos
    
    pair_names = cell(2,1);
    pair_names{1} = mammo_names{ii};
    pair_names{2} = mammo_names{ii};
    if strcmp(pair_names{1}(4), 'R')
        pair_names{2}(4) = 'L';
    else
        pair_names{2}(4) = 'R';
    end
    
    %Compute X-th percentile for maps of the abnormal and normal breast
    for jj = 1:2

        %load mask
        mask_name = dir([mask_dir '*' pair_names{jj}  '*.mat']);
        if length(mask_name) == 1
            mask = u_load([mask_dir mask_name(1).name]);
        else
            display(['Couldn''t load ' mask_dir '*' pair_names{jj} '*.mat']);
            continue;
        end
            
        for kk = 1:num_dists
            %load map
            map_name = dir([map_dir '*' pair_names{jj} '*' args.map_type '*' dist_names{kk} '*.mat']);
            if length(map_name) == 1
                map = load_uint8([map_dir map_name(1).name]);
            else
                display(['Couldn''t load ' map_dir '*' pair_names{jj} '*' dist_names{kk} '*.mat']);
                continue;
            end
            
            %smooth map be capture weighting
            if args.sigma
                map = imfilter(map, fspecial('gaussian', 5*args.sigma+1, args.sigma));
            end

            %check mask and map are same size
            if ~all(size(mask) == size(map))
                mask = imresize(mask, size(map));
            end

            %Work out the X-th percentile for the map
            percentiles(ii,jj,kk) = prctile(map(mask), args.percentile);
            clear map
        end
        clear mask
    end

end

%Workout threshhold as mean of the percentiles across all maps
thresh = squeeze(naNmean(naNmean(percentiles), 2));
