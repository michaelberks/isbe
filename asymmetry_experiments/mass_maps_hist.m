function [counts bins] = mass_maps_hist(data_type, map_dir, varargin)
%MASS_MAPS_ROC *Insert a one line summary here*
%   [] = mass_maps_maxima()
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
    'bin_width', 0.01,...
    'map_list', [],...
    'view', [],...
    'sigma', 0,...
    'mask_dir', 'masks',...
    'mammo_dir', 'mammograms',...
    'meta_dir', 'meta2',...
    'hit_method', 1,...
    'quiet', 1,...
    'num_mammos', inf,...
    'plot', 0);
clear varargin;

if args.quiet
    warning('off', 'load_uint8:missing_variables');
end

%form complete paths to map and mammogram directories
map_dir = [asymmetryroot 'data/' map_dir  '/' data_type '/'];
meta_dir = [asymmetryroot 'data/' args.mammo_dir  '/' data_type '/' args.meta_dir '/'];

if 0%exist(meta_dir, 'dir')
    abnormal = true;
else
    abnormal = false;
end

%get list of mammograms containing abnormalities
if isempty(args.map_list)
    map_list = dir([map_dir '/*' args.view '.mat']);
else
    map_list = args.map_list;
end
mammo_names = get_mammo_info(map_list);

if ~isempty(args.mask_dir);
    mask_dir = [asymmetryroot 'data/' args.mask_dir  '/' data_type '/'];
    [mask_names missing_idx] =...
        match_mammo_names(mask_dir, mammo_names);
    
    if ~isempty(missing_idx)
        display('warning: missing masks');
        map_list(missing_idx) = [];
        mammo_names(missing_idx) = [];
        mask_names(missing_idx) = [];
    end
end
    
num_mammos = min(args.num_mammos, length(mammo_names));
if abnormal
    %Match map names to mammogram meta data
    [meta_names missing_idx] =...
        match_mammo_names(meta_dir, mammo_names);

    %Check there's no missing maps
    if ~isempty(missing_idx)
        display('warning: missing maps');
        map_list(missing_idx) = [];
        meta_names(missing_idx) = [];
        num_mammos = min(args.num_mammos, length(map_list));
    end
end

%Initialise bins and counts - these will grow to match the range in each
%map
bins = 0;
counts = 0;

%For each map
for ii = 1:num_mammos
    
    %display user feedback
    display(['Processing image ' num2str(ii) ' of ' num2str(num_mammos)]);
    
    %load map
    map = load_uint8([map_dir map_list(ii).name]);
    [r c] = size(map);
    
    %smooth map
    if args.sigma
        map = imfilter(map, fspecial('gaussian', 5*args.sigma+1, args.sigma));
    end

    if args.plot && ii < 20
        figure; imagesc(map); axis image; colormap(jet(256)); hold on;
    end
    
    if ~isempty(args.mask_dir)
        map_mask = u_load([mask_dir mask_names{ii}]);
        map_mask = imresize(map_mask, [r c]);
    else
        map_mask = true(r, c);
    end
    
    if abnormal
        %load meta data
        meta = u_load([meta_dir meta_names{ii}]);
    
        %Update mask so only mass regions selected
        for jj = 1:length(meta)
        
            %Resize mass outlines
            meta_xy = meta{jj};
        
            %Take & of existing mask and mass region
            map_mask = map_mask & roipoly(map_mask, meta_xy(:,1)*c, meta_xy(:,2)*r);
        end
    end
    
    %Get min/max values of the map
    map_min = min(map(map_mask));   
    map_max = max(map(map_mask));
    
    %Extend the bins if necessary
    while map_min < bins(1)
        bins = [bins(1)-args.bin_width bins]; %#ok
        counts = [0 counts]; %#ok
    end
    while map_max > bins(end)
        bins = [bins bins(end)+args.bin_width]; %#ok
        counts = [counts 0]; %#ok
    end
    
    %Now get the histogram counts for this map
    map_counts = hist(map(map_mask), bins);
    
    %Add the set counts
    counts = counts + map_counts;
    
    %clear the map from memory
    clear map

end

