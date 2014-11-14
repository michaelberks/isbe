function [tp fp fp_pixels] = mass_maps_roc(map_dir, data_type, varargin)
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
    'thresh', linspace(0, 1, 101),...
    'map_list', [],...
    'map_type', [],...
    'sigma', 8,...
    'mask_dir', 'masks',...
    'mammo_dir', 'mammograms',...
    'meta_dir', 'meta2',...
    'exclusion_zone', 80,...
    'resize_factor', 0.5,...
    'hit_method', 3,...
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

if exist(meta_dir, 'dir')
    abnormal = true;
else
    abnormal = false;
end

%get list of mammograms containing abnormalities
if isempty(args.map_list)
    map_list = dir([map_dir '/*' args.map_type '*.mat']);
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
num_pts = length(args.thresh);

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
    
    %True positives can only occur for abnormal mammograms
    tp = zeros(num_mammos, num_pts);
else
    tp = [];
end

%pre-allocate false positives
fp = zeros(num_mammos, num_pts);
fp_pixels = zeros(num_mammos, num_pts);

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

    if abnormal
        %load meta data
        meta = u_load([meta_dir meta_names{ii}]);
    end
    
    if args.plot && ii < 20
        figure; imagesc(map); axis image; colormap(jet(256)); hold on;
    end
    
    if ~isempty(args.mask_dir)
        map_mask = load_uint8([mask_dir mask_names{ii}]);
        map_mask = imresize(map_mask, [r c]);
    else
        map_mask = true(r, c);
    end
    
    %Go through theshold points computing true/false positives
    for jj = 1:num_pts
        
        %theshold the map
        label_mask = (map >= args.thresh(jj)) & map_mask;
        
        %Get the connected components labels of the mask
        label = bwlabel(label_mask);
        num_components = max(label(:));

        if ~num_components
            break;
        end
        
        if abnormal
            %Go through each connected component and work it if its tp or a
            %false p
            hits = zeros(length(meta), 1);
            for kk = 1:length(meta)
        
                %Resize mass outlines
                meta_xy = meta{kk};
                
                if args.plot && jj == 1
                    plot(meta_xy(:,1)*c, meta_xy(:,2)*r, 'k');
                end
                
                %loop though each component
                for ll = 1:num_components
                    [label_y label_x] = find(label == ll);
                    hits(kk) = hits(kk) + ...
                        is_a_hit(label_x, label_y, meta_xy(:,1)*c, meta_xy(:,2)*r, map(label == ll), args.hit_method);
                end
            end
            %true positives are masses that are hits
            tp(ii,jj) = sum(hits > 0);
            
            %false positives are the number of components that weren't hits
            fp(ii, jj) = max(num_components - sum(hits), 0);
        else
            %fp count is number of components
            fp(ii, jj) = num_components;
            fp_pixels(ii, jj) = sum(label_mask(:));
        end
    end
    
    %clear the map from memory
    clear map

end

function [tf] = is_a_hit(label_x, label_y, mass_x, mass_y, label_vals, method)

switch method
    case 1
        %Intersection of regions is greater than 50% of label size
        tf = sum(inpolygon(label_x, label_y, mass_x, mass_y)) > 0.5*length(label_x(:));
        
    case 2
        %Intersection of regions is greater than 50% of label size & mass
        %size
        tf = sum(inpolygon(label_x, label_y, mass_x, mass_y)) > 0.5*length(label_x(:)) &&...
             sum(inpolygon(label_x, label_y, mass_x, mass_y)) > 0.5*length(mass_x(:)) ;
    
    case 3
        %Who cares, as long as there's some overlap...
        tf = any(inpolygon(label_x, label_y, mass_x, mass_y));
        
    case 4
        %Centre-of-mass of label (with label values as weights) must lie
        %inside mass
        label_sum = sum(label_vals);
        x_mean = sum(label_x .* label_vals) / label_sum;
        y_mean = sum(label_y .* label_vals) / label_sum;
        tf = inpolygon(x_mean, y_mean, mass_x, mass_y);
        
    case 5
        %Centre-of-mass of label (with uniform weights) must lie
        %inside mass
        tf = inpolygon(mean(label_x), mean(label_y), mass_x, mass_y);
end

