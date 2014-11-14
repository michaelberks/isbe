function [hits num_maxima maxima_ranks maxima_scores] = mass_maps_maxima(map_dir, thresh, varargin)
%MASS_MAPS_MAXIMA *Insert a one line summary here*
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
    'mammo_dir', 'mammograms',...
    'mammo_list', [],...
    'data_type', '2004_screening/abnormals',...
    'sigma', 0,...
    'mask_dir', 'masks',...
    'meta_dir', 'meta2',...
    'map_type', [],...
    'dist_ranges', [],...
    'exclusion_zone', 80,...
    'quiet', 1,...
    'plot', 0);
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

hits = nan(num_mammos,2,num_dists);
num_maxima = nan(num_mammos,2,num_dists);
maxima_scores = nan(num_mammos,2,num_dists);
maxima_ranks = nan(num_mammos,2,num_dists);

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
    
    if args.plot
        figs = zeros(num_dists);
        for jj = 1:num_dists
            figs(jj) = figure('Name', [pair_names{1}(1:3) ' ' pair_names{1}(5:6) ' Distance weighting = ' num2str(args.dist_ranges(jj))]);
        end
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
        
        if jj == 1
            %load mass outline for the abnormal mammo
            meta_name = dir([mammo_dir '/' args.meta_dir '/*' pair_names{jj}  '*.mat']);
            if length(mask_name) == 1
                meta = u_load([mammo_dir '/' args.meta_dir '/' meta_name(1).name]);
            else
                display(['Couldn''t load ' mask_dir '*' pair_names{jj} '*.mat']);
                continue;
            end
        else
            %flip the mask for the normal mammogram
            mask = fliplr(mask);
        end
        
        for kk = 1:num_dists
            %load map
            map_name = dir([map_dir '*' pair_names{jj} '*' args.map_type '*' dist_names{kk} '*.mat']);
            if length(map_name) == 1
                map = load_uint8([map_dir map_name(1).name]);
                [r c] = size(map);
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
                mask = imresize(mask, [r c]);
            end
            
            if jj == 2
                %flip the normal map
                map = fliplr(map);
            end
            
            %get positions and values of local maxima in map
            [maxima_pos maxima_vals] = local_image_maxima(map, args.exclusion_zone, mask, thresh(kk));
            
            %Count number maxima
            num_maxima(ii,jj,kk) = size(maxima_pos,1);

            %Work out whether the maxima lie within the mass outline
            map_hits = [];
            for ll = 1:length(meta)
                meta_xy = meta{ll};
                map_hits = [map_hits; inpolygon(maxima_pos(:,1), maxima_pos(:,2), meta_xy(:,1)*c, meta_xy(:,2)*r)]; %#ok
            end

            %if we have any hits, get their rank and value
            if any(map_hits)
                hits(ii,jj,kk) = 1;
                maxima_ranks(ii,jj,kk) = rem(find(map_hits, 1)-1, length(maxima_vals))+1;
                maxima_scores(ii,jj,kk) = maxima_vals(maxima_ranks(ii,jj,kk));
            else
                hits(ii,jj,kk) = 0;
            end
            
            if args.plot
                figure(figs(kk));
                subplot(1,2,jj); imagesc(map); axis image; colormap(jet(256)); hold on;
                plot(maxima_pos(:,1), maxima_pos(:,2), 'kx');
                for ll = 1:length(meta)
                    plot(meta{ll}(:,1)*c, meta{ll}(:,2)*r, 'k');
                end
                if hits(ii,jj,kk)
                    title('It''s a hit!');
                else
                    title('It''s a miss!');
                end
                if jj == 2
                    subplot(1,2,1);
                    mm = caxis;
                    mm(2) = max(mm(2), max(map(:)));
                    caxis(mm);
                    subplot(1,2,2);
                    caxis(mm);
                end
                    
            end
            
            
            clear map
        end
        clear mask
    end

end