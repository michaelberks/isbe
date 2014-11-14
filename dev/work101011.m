rad064_list = dir('C:\isbe\asymmetry_project\data\radial_maps\2004_screening\abnormals\*064.mat');
rad_names = get_mammo_info(rad064_list);
pcn064 = zeros(length(rad064_list),5);
pts = [10 25 50 75 90];
for ii = 1:length(rad064_list)
    rad_map = u_load(['C:\isbe\asymmetry_project\data\radial_maps\2004_screening\abnormals\' rad064_list(ii).name]);
    mask = u_load(['C:\isbe\asymmetry_project\data\masks\2004_screening\abnormals\o04_' rad_names{ii} '_mask.mat']);
    mask = imresize(mask, size(rad_map));
    pcn064(ii,:) = prctile(rad_map(mask), pts);
    clear rad_map masks
end
%
rad128_list = dir('C:\isbe\asymmetry_project\data\radial_maps\2004_screening\abnormals\*128.mat');
rad_names = get_mammo_info(rad128_list);
pcn128 = zeros(length(rad128_list),5);
pts = [10 25 50 75 90];
for ii = 1:length(rad128_list)
    rad_map = u_load(['C:\isbe\asymmetry_project\data\radial_maps\2004_screening\abnormals\' rad128_list(ii).name]);
    mask = u_load(['C:\isbe\asymmetry_project\data\masks\2004_screening\abnormals\o04_' rad_names{ii} '_mask.mat']);
    mask = imresize(mask, size(rad_map));
    pcn128(ii,:) = prctile(rad_map(mask), pts);
    clear rad_map masks
end
rad256_list = dir('C:\isbe\asymmetry_project\data\radial_maps\2004_screening\abnormals\*256.mat');
rad_names = get_mammo_info(rad256_list);
pcn256 = zeros(length(rad256_list),5);
pts = [10 25 50 75 90];
for ii = 1:length(rad256_list)
    rad_map = u_load(['C:\isbe\asymmetry_project\data\radial_maps\2004_screening\abnormals\' rad256_list(ii).name]);
    mask = u_load(['C:\isbe\asymmetry_project\data\masks\2004_screening\abnormals\o04_' rad_names{ii} '_mask.mat']);
    mask = imresize(mask, size(rad_map));
    pcn256(ii,:) = prctile(rad_map(mask), pts);
    clear rad_map masks
end
%%
hits = zeros(292,3);

for jj = 1:3
    rad_list = dir(['C:\isbe\asymmetry_project\data\radial_maps\2004_screening\abnormals\*' zerostr(2^(5+jj),3) '.mat']);
    rad_names = get_mammo_info(rad_list);


    for ii = 1:length(rad_list)
        if ~exist(['C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\mat\meta\' rad_names{ii} '_meta.mat'], 'file')
            continue;
        end

        rad_map = u_load(['C:\isbe\asymmetry_project\data\radial_maps\2004_screening\abnormals\' rad_list(ii).name]);
        rad_map = imfilter(rad_map, fspecial('gaussian', 40, 8));
        mask = u_load(['C:\isbe\asymmetry_project\data\masks\2004_screening\abnormals\o04_' rad_names{ii} '_mask.mat']);
        mask = imresize(mask, size(rad_map));
        meta = u_load(['C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\mat\meta\' rad_names{ii} '_meta.mat']);

        meta = meta / 2;
        y = prctile(rad_map(mask), 90);
        [maxima_pos] = local_image_maxima(rad_map, 80, mask, y);
        hits(ii,jj) = any(inpolygon(maxima_pos(:,1), maxima_pos(:,2), meta(:,1), meta(:,2)));

    %     figure; imagesc(rad_map); axis image; colormap(jet(256)); hold on;
    %     plot(maxima_pos(:,1), maxima_pos(:,2), 'ko');
    %     plot(meta(:,1), meta(:,2), 'k');
    %     
    %     if hits(ii)
    %         title('It''s a score!');
    %     else
    %         title('It''s a miss!');
    %     end

        clear rad_map masks meta
    end
end   
%%
hits_u = zeros(292,5);

for jj = 1:5
    rad_list = dir(['C:\isbe\asymmetry_project\data\radial_maps\2004_screening\abnormals\*' zerostr(2^(3+jj),3) '.mat']);
    rad_names = get_mammo_info(rad_list);


    for ii = 1:length(rad_list)
        if ~exist(['C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\mat\meta\' rad_names{ii} '_meta.mat'], 'file')
            continue;
        end

        rad_map = u_load(['C:\isbe\asymmetry_project\data\radial_maps\2004_screening\abnormals\' rad_list(ii).name]);
        rad_map = imfilter(rad_map, fspecial('gaussian', 40, 8));
        mask = u_load(['C:\isbe\asymmetry_project\data\masks\2004_screening\abnormals\o04_' rad_names{ii} '_mask.mat']);
        mask = imresize(mask, size(rad_map));
        meta = u_load(['C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\mat\meta\' rad_names{ii} '_meta.mat']);

        meta = meta / 2;
        y = prctile(rad_map(mask), 90);
        [maxima_pos] = local_image_maxima(rad_map, 80, mask, y);
        hits_u(ii,jj) = any(inpolygon(maxima_pos(:,1), maxima_pos(:,2), meta(:,1), meta(:,2)));

%         figure; imagesc(rad_map); axis image; colormap(jet(256)); hold on;
%         plot(maxima_pos(:,1), maxima_pos(:,2), 'ko');
%         plot(meta(:,1), meta(:,2), 'k');
%         
%         if hits_u(ii)
%             title('It''s a score!');
%         else
%             title('It''s a miss!');
%         end

        clear rad_map masks meta
    end
end    
%%
hits_w = zeros(292,5);

for jj = 3:5
    rad_list = dir(['C:\isbe\asymmetry_project\data\weighted_radial_maps\2004_screening\abnormals\*' zerostr(2^(3+jj),3) '.mat']);
    rad_names = get_mammo_info(rad_list);


    for ii = 1:length(rad_list)
        if ~exist(['C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\mat\meta\' rad_names{ii} '_meta.mat'], 'file')
            continue;
        end

        rad_map = u_load(['C:\isbe\asymmetry_project\data\weighted_radial_maps\2004_screening\abnormals\' rad_list(ii).name]);
        rad_map = imfilter(rad_map, fspecial('gaussian', 40, 8));
        mask = u_load(['C:\isbe\asymmetry_project\data\masks\2004_screening\abnormals\o04_' rad_names{ii} '_mask.mat']);
        mask = imresize(mask, size(rad_map));
        meta = u_load(['C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\mat\meta\' rad_names{ii} '_meta.mat']);

        meta = meta / 2;
        y = prctile(rad_map(mask), 90);
        [maxima_pos] = local_image_maxima(rad_map, 80, mask, y);
        hits_w(ii,jj) = any(inpolygon(maxima_pos(:,1), maxima_pos(:,2), meta(:,1), meta(:,2)));

%         figure; imagesc(rad_map); axis image; colormap(jet(256)); hold on;
%         plot(maxima_pos(:,1), maxima_pos(:,2), 'ko');
%         plot(meta(:,1), meta(:,2), 'k');
%         
%         if hits_w(ii)
%             title('It''s a score!');
%         else
%             title('It''s a miss!');
%         end

        clear rad_map masks meta
    end
end   
%%
hits_nw = zeros(292,5);

for jj = 3:5
    rad_list = dir(['C:\isbe\asymmetry_project\data\weighted_radial_maps\2004_screening\abnormals\*' zerostr(2^(3+jj),3) '.mat']);
    rad_names = get_mammo_info(rad_list);


    for ii = 1:length(rad_list)
        
        if ~exist(['C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\mat\meta\' rad_names{ii} '_meta.mat'], 'file')
            continue;
        end
        
        mammo_name = rad_names{ii};
        rad_name = rad_list(ii).name;
        if strcmp(mammo_name(4), 'R')
            mammo_name(4) = 'L';
            rad_name(4) = 'L';
        else
            mammo_name(4) = 'R';
            rad_name(4) = 'R';
        end
        
        rad_map = u_load(['C:\isbe\asymmetry_project\data\weighted_radial_maps\2004_screening\abnormals\' rad_name]);
        rad_map = fliplr(imfilter(rad_map, fspecial('gaussian', 40, 8)));
        mask = u_load(['C:\isbe\asymmetry_project\data\masks\2004_screening\abnormals\o04_' mammo_name '_mask.mat']);
        mask = fliplr(imresize(mask, size(rad_map)));
        meta = u_load(['C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\mat\meta\' rad_names{ii} '_meta.mat']);

        meta = meta / 2;
        y = prctile(rad_map(mask), 90);
        [maxima_pos] = local_image_maxima(rad_map, 80, mask, y);
        hits_nw(ii,jj) = any(inpolygon(maxima_pos(:,1), maxima_pos(:,2), meta(:,1), meta(:,2)));

%         figure; imagesc(rad_map); axis image; colormap(jet(256)); hold on;
%         plot(maxima_pos(:,1), maxima_pos(:,2), 'ko');
%         plot(meta(:,1), meta(:,2), 'k');
%         
%         if hits_nw(ii,jj)
%             title('It''s a score!');
%         else
%             title('It''s a miss!');
%         end

        clear rad_map masks meta
    end
end  
%%
%**************************************************************************
%**************************************************************************
mammo_list = dir('C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\mat\meta\*.mat');
mammo_names = get_mammo_info(mammo_list);

map_dir = 'weighted_radial_maps';

percentile_95_n = nan(146,6);
percentile_95_a = nan(146,6);

for ii = 1:146
    
    contra_name = mammo_names{ii};

    if strcmp(contra_name(4), 'R')
        contra_name(4) = 'L';
    else
        contra_name(4) = 'R';
    end
    
    for jj = 1:6
        
        %Compute 95th percentile for the abnormal breast
        try
            rad_map = u_load(['C:\isbe\asymmetry_project\data\' map_dir '\2004_screening\abnormals\' mammo_names{ii} '_rad_map_' zerostr(2^(2+jj),3) '.mat']);
        catch
            display(['Couldn''t load C:\isbe\asymmetry_project\data\' map_dir '\2004_screening\abnormals\' mammo_names{ii} '_rad_map_' zerostr(2^(2+jj),3) '.mat']);
            continue;
        end

        mask = u_load(['C:\isbe\asymmetry_project\data\masks\2004_screening\abnormals\o04_' mammo_names{ii} '_mask.mat']);
        mask = imresize(mask, size(rad_map));
        
        percentile_95_a(ii,jj) = prctile(rad_map(mask), 95);
        clear rad_map mask
        
        %Compute 95th percentile for the normal breast
        try
            rad_map = u_load(['C:\isbe\asymmetry_project\data\' map_dir '\2004_screening\abnormals\' contra_name '_rad_map_' zerostr(2^(2+jj),3) '.mat']);
        catch
            display(['Couldn''t load C:\isbe\asymmetry_project\data\' map_dir '\2004_screening\abnormals\' contra_name '_rad_map_' zerostr(2^(2+jj),3) '.mat']);
            continue;
        end

        mask = u_load(['C:\isbe\asymmetry_project\data\masks\2004_screening\abnormals\o04_' contra_name '_mask.mat']);
        mask = imresize(mask, size(rad_map));
        
        percentile_95_n(ii,jj) = prctile(rad_map(mask), 95);
        clear rad_map mask
    end
end

thresh_a = naNmean(percentile_95_a); display(thresh_a);
thresh_n = naNmean(percentile_95_n); display(thresh_n);
thresh = naNmean([percentile_95_a; percentile_95_n]);

save(['C:\isbe\asymmetry_project\experiments\radial_maps\thresh_' map_dir '.mat'], 'percentile_95*', 'thresh*');
%%

map_dir = 'weighted_radial_maps';
load(['C:\isbe\asymmetry_project\experiments\radial_maps\thresh_' map_dir '.mat'], 'thresh*');

mammo_list = dir('C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\mat\meta\*.mat');
mammo_names = get_mammo_info(mammo_list);

hits_n = nan(146,6);
hits_a = nan(146,6);

maxima_n = nan(146,6);
maxima_a = nan(146,6);

score_n = nan(146,6);
score_a = nan(146,6);

rank_n = nan(146,6);
rank_a = nan(146,6);
%
for ii = 1:146
    
    contra_name = mammo_names{ii};

    if strcmp(contra_name(4), 'R')
        contra_name(4) = 'L';
    else
        contra_name(4) = 'R';
    end
    
    for jj = 1:6
        
        %Compute hits for the abnormal breast
        try %and load radial map
            rad_map = u_load(['C:\isbe\asymmetry_project\data\' map_dir '\2004_screening\abnormals\' mammo_names{ii} '_rad_map_' zerostr(2^(2+jj),3) '.mat']);
        catch
            display(['Couldn''t load C:\isbe\asymmetry_project\data\' map_dir '\2004_screening\abnormals\' mammo_names{ii} '_rad_map_' zerostr(2^(2+jj),3) '.mat']);
            continue;
        end

        %Load mask of breast and resize
        mask = u_load(['C:\isbe\asymmetry_project\data\masks\2004_screening\abnormals\o04_' mammo_names{ii} '_mask.mat']);
        mask = imresize(mask, size(rad_map));
        
        %Load mass outline in breast and resize
        meta = u_load(['C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\mat\meta2\' mammo_names{ii} '_meta.mat']);
        

        %Workout if we have a hit
        [maxima_pos maxima_vals] = local_image_maxima(rad_map, 80, mask, thresh(jj));
        maxima_a(ii,jj) = size(maxima_pos,1); %Count number maxima
        
        hits = [];
        for kk = 1:length(meta)
            meta{kk} = meta{kk} / 2;
            hits = [hits; inpolygon(maxima_pos(:,1), maxima_pos(:,2), meta{kk}(:,1), meta{kk}(:,2))]; %#ok
        end
        
        if any(hits)
            hits_a(ii,jj) = 1;
            rank_a(ii,jj) = rem(find(hits, 1)-1, length(maxima_vals))+1;
            score_a(ii,jj) = maxima_vals(rank_a(ii,jj));
        end    
        clear rad_map mask
        
         %Compute hits for the normal breast
        try %and load radial map
            rad_map = fliplr(u_load(['C:\isbe\asymmetry_project\data\' map_dir '\2004_screening\abnormals\' contra_name '_rad_map_' zerostr(2^(2+jj),3) '.mat']));
        catch
            display(['Couldn''t load C:\isbe\asymmetry_project\data\' map_dir '\2004_screening\abnormals\' contra_name '_rad_map_' zerostr(2^(2+jj),3) '.mat']);
            continue;
        end

        %Load mask of breast and resize
        mask = u_load(['C:\isbe\asymmetry_project\data\masks\2004_screening\abnormals\o04_' contra_name '_mask.mat']);
        mask = fliplr(imresize(mask, size(rad_map)));

        %Workout if we have a hit
        [maxima_pos maxima_vals] = local_image_maxima(rad_map, 80, mask, thresh(jj));
        maxima_n(ii,jj) = size(maxima_pos,1);
        
        hits = [];
        for kk = 1:length(meta)
            hits = [hits; inpolygon(maxima_pos(:,1), maxima_pos(:,2), meta{kk}(:,1), meta{kk}(:,2))]; %#ok
        end
        
        if any(hits)
            hits_n(ii,jj) = 1;
            rank_n(ii,jj) = rem(find(hits, 1)-1, length(maxima_vals))+1;
            score_n(ii,jj) = maxima_vals(rank_n(ii,jj));
        end    
        clear rad_map mask
    end
end
save(['C:\isbe\asymmetry_project\experiments\radial_maps\hits_' map_dir '.mat'], 'hits*', 'score*', 'maxima*', 'rank*');
%%
[sum(nansum(hits_a,2)==0) sum(nansum(hits_a,2)==0)]
[sum(nansum(hits_a,2)==1) sum(nansum(hits_n,2)==1)]
[sum(nansum(hits_a,2)==2) sum(nansum(hits_n,2)==2)]
[sum(nansum(hits_a,2)==3) sum(nansum(hits_n,2)==3)]
[sum(nansum(hits_a,2)==4) sum(nansum(hits_n,2)==4)]
[sum(nansum(hits_a,2)==5) sum(nansum(hits_n,2)==5)]
[sum(nansum(hits_a,2)==6) sum(nansum(hits_n,2)==6)]

 
%%
rad_list = [...
    dir('C:\isbe\asymmetry_project\data\radial_maps\2004_screening\abnormals\*008.mat');...
    dir('C:\isbe\asymmetry_project\data\radial_maps\2004_screening\abnormals\*016.mat');...
    dir('C:\isbe\asymmetry_project\data\radial_maps\2004_screening\abnormals\*032.mat')];
for ii = 1:length(rad_list)
    copyfile(...
        ['C:\isbe\asymmetry_project\data\radial_maps\2004_screening\abnormals\' rad_list(ii).name],...
        ['C:\isbe\asymmetry_project\data\weighted_radial_maps\2004_screening\abnormals\' rad_list(ii).name]);
end
%%
mammo_list = dir('C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\mat\meta\*.mat');
mammo_names = get_mammo_info(mammo_list);

map_dir = 'mass_maps';

percentile_95_n = nan(146,1);
percentile_95_a = nan(146,1);
%
for ii = 1:146
    
    contra_name = mammo_names{ii};

    if strcmp(contra_name(4), 'R')
        contra_name(4) = 'L';
    else
        contra_name(4) = 'R';
    end
        
    %Compute 95th percentile for the abnormal breast
    try
        mass_map = u_load(['C:\isbe\asymmetry_project\data\' map_dir '\2004_screening\abnormals\o04_' mammo_names{ii} '_mass.mat']);
    catch
        display(['Couldn''t load C:\isbe\asymmetry_project\data\' map_dir '\2004_screening\abnormals\o04_' mammo_names{ii} '_mass.mat']);
        continue;
    end

    mask = u_load(['C:\isbe\asymmetry_project\data\masks\2004_screening\abnormals\o04_' mammo_names{ii} '_mask.mat']);
    mask = imresize(mask, size(mass_map));

    percentile_95_a(ii,1) = prctile(mass_map(mask), 95);
    clear mass_map mask

    %Compute 95th percentile for the normal breast
    try
        mass_map = u_load(['C:\isbe\asymmetry_project\data\' map_dir '\2004_screening\abnormals\o04_' contra_name '_mass.mat']);
    catch
        display(['Couldn''t load C:\isbe\asymmetry_project\data\' map_dir '\2004_screening\abnormals\o04_' contra_name '_mass.mat']);
        continue;
    end

    mask = u_load(['C:\isbe\asymmetry_project\data\masks\2004_screening\abnormals\o04_' contra_name '_mask.mat']);
    mask = imresize(mask, size(mass_map));

    percentile_95_n(ii,1) = prctile(mass_map(mask), 95);
    clear mass_map mask

end

thresh_a = naNmean(percentile_95_a); display(thresh_a);
thresh_n = naNmean(percentile_95_n); display(thresh_n);
thresh = naNmean([percentile_95_a; percentile_95_n]);

save(['C:\isbe\asymmetry_project\experiments\radial_maps\thresh_' map_dir '.mat'], 'percentile_95*', 'thresh*');
%
map_dir = 'mass_maps';

mammo_list = dir('C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\mat\meta\*.mat');
mammo_names = get_mammo_info(mammo_list);

hits_n = nan(146,1);
hits_a = nan(146,1);

maxima_n = nan(146,1);
maxima_a = nan(146,1);

score_n = nan(146,1);
score_a = nan(146,1);

rank_n = nan(146,1);
rank_a = nan(146,1);
%
do_plot = 0;

for ii = 1:146
    %ii = idx(iii);
    contra_name = mammo_names{ii};

    if strcmp(contra_name(4), 'R')
        contra_name(4) = 'L';
    else
        contra_name(4) = 'R';
    end
        
    %Compute hits for the abnormal breast
    try %and load radial map
        mass_map = u_load(['C:\isbe\asymmetry_project\data\' map_dir '\2004_screening\abnormals\o04_' mammo_names{ii} '_mass.mat']);
    catch
        display(['Couldn''t load C:\isbe\asymmetry_project\data\' map_dir '\2004_screening\abnormals\o04_' mammo_names{ii} '_mass.mat']);
        continue;
    end

    %Load mask of breast and resize
    mask = u_load(['C:\isbe\asymmetry_project\data\masks\2004_screening\abnormals\o04_' mammo_names{ii} '_mask.mat']);
    mask = imresize(mask, size(mass_map));

    %Load mass outline in breast and resize
    meta = u_load(['C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\mat\meta2\' mammo_names{ii} '_meta.mat']);


    %Workout if we have a hit
    [maxima_pos maxima_vals] = local_image_maxima(mass_map, 80, mask, thresh);
    maxima_a(ii,1) = size(maxima_pos,1); %Count number maxima

    if do_plot
        figure; subplot(1,2,1); imagesc(mass_map); axis image; colormap(jet(256)); hold on;
        plot(maxima_pos(:,1), maxima_pos(:,2), 'ko');
    end
    
    hits = [];
    for kk = 1:length(meta)
        meta{kk} = meta{kk} / 2;
        hits = [hits; inpolygon(maxima_pos(:,1), maxima_pos(:,2), meta{kk}(:,1), meta{kk}(:,2))]; %#ok
        
        if do_plot
            plot(meta{kk}(:,1), meta{kk}(:,2), 'k');
        end
    end

    if any(hits)
        hits_a(ii,1) = 1;
        rank_a(ii,1) = rem(find(hits, 1)-1, length(maxima_vals))+1;
        score_a(ii,1) = maxima_vals(rank_a(ii,1));
    end    
    clear mass_map mask

     %Compute hits for the normal breast
    try %and load radial map
        mass_map = fliplr(u_load(['C:\isbe\asymmetry_project\data\' map_dir '\2004_screening\abnormals\o04_' contra_name '_mass.mat']));
    catch
        display(['Couldn''t load C:\isbe\asymmetry_project\data\' map_dir '\2004_screening\abnormals\o04_' contra_name '_mass.mat']);
        continue;
    end

    %Load mask of breast and resize
    mask = u_load(['C:\isbe\asymmetry_project\data\masks\2004_screening\abnormals\o04_' contra_name '_mask.mat']);
    mask = fliplr(imresize(mask, size(mass_map)));

    %Workout if we have a hit
    [maxima_pos maxima_vals] = local_image_maxima(mass_map, 80, mask, thresh);
    maxima_n(ii,1) = size(maxima_pos,1);

    if do_plot
        subplot(1,2,2); imagesc(mass_map); axis image; colormap(jet(256)); hold on;
        plot(maxima_pos(:,1), maxima_pos(:,2), 'ko');
    end
    
    hits = [];
    for kk = 1:length(meta)
        hits = [hits; inpolygon(maxima_pos(:,1), maxima_pos(:,2), meta{kk}(:,1), meta{kk}(:,2))]; %#ok
        if do_plot
            plot(meta{kk}(:,1), meta{kk}(:,2), 'k');
        end
    end

    if any(hits)
        hits_n(ii,1) = 1;
        rank_n(ii,1) = rem(find(hits, 1)-1, length(maxima_vals))+1;
        score_n(ii,1) = maxima_vals(rank_n(ii,1));
    end    
    clear mass_map mask
    
end
save(['C:\isbe\asymmetry_project\experiments\radial_maps\hits_' map_dir '.mat'], 'hits*', 'score*', 'maxima*', 'rank*');
%%
meta_list = dir('C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\mat\meta2\*.mat');
areas = 0;
for ii = length(meta_list)
    load(['C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\mat\meta2\' meta_list(ii).name]);
    
    for jj = 1:length(meta_xy)
        areas = areas + polyarea(meta_xy{jj}(:,1), meta_xy{jj}(:,2)) / 4;
    end
end
%%
[training_data training_labels] = sample_mass_training_data(... % non-strict mode
    'num_samples', 2e3,...
    'abnormal_data', '2004_screening/abnormals/',...
    'normal_data', '2004_screening/normals/',... 
    'image_dir', [asymmetryroot, 'data/mammograms/'],...
    'mask_dir', [asymmetryroot, 'data/mass_masks/'],...
    'radial_dir', [asymmetryroot, 'data/weighted_radial_maps/'],...
    'mass_dir', [asymmetryroot, 'data/mass_maps/'],...
    'fold_id', 1,...
    'num_folds', 10,...
    'view', [],...
    'dist_range', [32 64 128 256],...
    'sigma_range', [1 2 4 8],...
    'angular_res', 1,...
    'save_path', []);
%%
lines1 = dir('C:\isbe\asymmetry_project\data\line_maps\2004_screening\abnormals\*.mat');
lines2 = dir('Z:\asymmetry_project\data\mammograms\2004_screening\abnormals\mat\results\191905\*.mat');

for ii = 4
    map1 = u_load(['C:\isbe\asymmetry_project\data\line_maps\2004_screening\abnormals\' lines1(ii).name]);
    map2 = u_load(['Z:\asymmetry_project\data\mammograms\2004_screening\abnormals\mat\results\191905\' lines2(ii).name]);
    
    figure;
    a1 = subplot(1,2,1); imagesc(map1>.5); axis image;
    a2 = subplot(1,2,2); imagesc(map2>.5); axis image;
    linkaxes([a1 a2]);
    clear map*
end
%%
lines1 = dir('C:\isbe\asymmetry_project\data\orientation_maps\2004_screening\abnormals\*.mat');
lines2 = dir('Z:\asymmetry_project\data\mammograms\2004_screening\abnormals\mat\results\191934\*.mat');

for ii = 1:length(lines1)
    copyfile(...
        ['Z:\asymmetry_project\data\mammograms\2004_screening\abnormals\mat\results\191934\' lines2(ii).name],...
        ['C:\isbe\asymmetry_project\data\orientation_maps\2004_screening\abnormals\' lines1(ii).name]);
    copyfile(...
        ['Z:\asymmetry_project\data\mammograms\2004_screening\abnormals\mat\results\191934\' lines2(ii).name],...
        ['Z:\asymmetry_project\data\orientation_maps\2004_screening\abnormals\' lines1(ii).name]);
end
%%
lines1 = dir('C:\isbe\asymmetry_project\data\orientation_maps\2004_screening\normals\*.mat');
for ii = 1:length(lines1)
    p_name = ['Z:\asymmetry_project\data\mammograms\2004_screening\normals\mat\results\191934\probability_image' zerostr(ii,3) '.mat'];
    if exist(p_name, 'file')
        movefile(...
            p_name,...
            ['Z:\asymmetry_project\data\mammograms\2004_screening\normals\mat\results\191934\' lines1(ii).name]);
    end

end
%%
clear
lines1 = dir('C:\isbe\asymmetry_project\data\line_maps\2004_screening\normals\*.mat');
lines2 = dir('Z:\asymmetry_project\data\mammograms\2004_screening\normals\mat\results\191905\*.mat');

for ii = 1:length(lines1)
    copyfile(...
        ['Z:\asymmetry_project\data\mammograms\2004_screening\normals\mat\results\191905\' lines2(ii).name],...
        ['C:\isbe\asymmetry_project\data\line_maps\2004_screening\normals\' lines1(ii).name]);
    copyfile(...
        ['Z:\asymmetry_project\data\mammograms\2004_screening\normals\mat\results\191905\' lines2(ii).name],...
        ['Z:\asymmetry_project\data\line_maps\2004_screening\normals\' lines1(ii).name]);
end
%%
clear
lines1 = dir('C:\isbe\asymmetry_project\data\line_maps\2004_screening\normals\*.mat');
lines2 = dir('Z:\asymmetry_project\data\mammograms\2004_screening\normals\mat\results\191934\*.mat');

for ii = 1:length(lines1)
    copyfile(...
        ['Z:\asymmetry_project\data\mammograms\2004_screening\normals\mat\results\191934\' lines2(ii).name],...
        ['C:\isbe\asymmetry_project\data\orientation_maps\2004_screening\normals\' lines1(ii).name]);
    copyfile(...
        ['Z:\asymmetry_project\data\mammograms\2004_screening\normals\mat\results\191934\' lines2(ii).name],...
        ['Z:\asymmetry_project\data\orientation_maps\2004_screening\normals\' lines1(ii).name]);
end