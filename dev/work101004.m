ab_sum = 0;
ab_N = 0;
norm_sum = 0;
norm_N = 0;
ab_means = [];
norm_means = [];

for fold = 1:10
    
    ab_list = dir(['Z:\asymmetry_project\data\mammograms\2004_screening\contralateral_roi\abnormals\results\roi_' zerostr(fold, 2) '\*.mat']);
    norm_list = dir(['Z:\asymmetry_project\data\mammograms\2004_screening\contralateral_roi\normals\results\roi_' zerostr(fold, 2) '\*.mat']);
    
    for ii = 1:length(norm_list)
        
        ab_map = u_load(['Z:\asymmetry_project\data\mammograms\2004_screening\contralateral_roi\abnormals\results\roi_' zerostr(fold, 2) '\' ab_list(ii).name]);
        norm_map = u_load(['Z:\asymmetry_project\data\mammograms\2004_screening\contralateral_roi\normals\results\roi_' zerostr(fold, 2) '\' norm_list(ii).name]);
        
        if strcmpi(ab_list(ii).name(4), 'R')
            ab_map = fliplr(ab_map);
        else
            norm_map = fliplr(norm_map);
        end
        ab_mask = u_load(['Z:\asymmetry_project\data\masks\2004_screening\contralateral_roi\abnormals\' ab_list(ii).name(1:6) '_mask.mat']);
        norm_mask = u_load(['Z:\asymmetry_project\data\masks\2004_screening\contralateral_roi\normals\' norm_list(ii).name(1:6) '_mask.mat']);       
        
        ab_im_sum = sum(1 - ab_map(ab_mask));
        ab_im_N = sum(ab_mask(:));
        norm_im_sum = sum(1 - norm_map(norm_mask));
        norm_im_N = sum(norm_mask(:));
        
        %display(['Image ' num2str(ii) ': ' num2str([ab_im_sum/ab_im_N norm_im_sum/norm_im_N], 3) ]); 
        
        ab_sum = ab_sum + ab_im_sum;
        ab_N = ab_N + ab_im_N;
        norm_sum = norm_sum + norm_im_sum;
        norm_N = norm_N + norm_im_N;
        
        ab_means(end+1) = ab_im_sum/ab_im_N; %#ok
        norm_means(end+1) = norm_im_sum/norm_im_N; %#ok
        
    end
end
display(ab_sum / ab_N);
display(norm_sum / norm_N);
[h p] = ttest(ab_means, norm_means) %#ok
sum(ab_means > norm_means) %#ok
%%
fold = 10;
    
ab_list = dir(['Z:\asymmetry_project\data\mammograms\2004_screening\contralateral_roi\abnormals\results\roi_' zerostr(fold, 2) '\*.mat']);
norm_list = dir(['Z:\asymmetry_project\data\mammograms\2004_screening\contralateral_roi\normals\results\roi_' zerostr(fold, 2) '\*.mat']);

for ii = 1:length(ab_list)

    ab_map = u_load(['Z:\asymmetry_project\data\mammograms\2004_screening\contralateral_roi\abnormals\results\roi_' zerostr(fold, 2) '\' ab_list(ii).name]);
    norm_map = u_load(['Z:\asymmetry_project\data\mammograms\2004_screening\contralateral_roi\normals\results\roi_' zerostr(fold, 2) '\' norm_list(ii).name]);
    
    ab_image = u_load(['Z:\asymmetry_project\data\mammograms\2004_screening\contralateral_roi\abnormals\' ab_list(ii).name(1:6)]);
    norm_image = u_load(['Z:\asymmetry_project\data\mammograms\2004_screening\contralateral_roi\normals\' norm_list(ii).name(1:6)]);

    min_prob = min(min(ab_map(:)), min(norm_map(:)));
    max_prob = max(max(ab_map(:)), max(norm_map(:)));
    
    if strcmpi(ab_list(ii).name(4), 'R')
        ab_map = fliplr(ab_map);
    else
        norm_map = fliplr(norm_map);
    end
    
    figure; 
    subplot(2,2,1); imagesc(ab_image); axis image; colormap(jet(256));
    subplot(2,2,2); imagesc(norm_image); axis image;
    subplot(2,2,3); imagesc(1-ab_map); axis image; caxis([0 1]);
    subplot(2,2,4); imagesc(1-norm_map); axis image; caxis([0 1]);

end
%%    
fold = 4;
    
ab_list = dir(['Z:\asymmetry_project\data\mammograms\2004_screening\contralateral\abnormals\results\roi_' zerostr(fold, 2) '\*.mat']);
norm_list = dir(['Z:\asymmetry_project\data\mammograms\2004_screening\contralateral\normals\results\roi_' zerostr(fold, 2) '\*.mat']);

for ii = 1:length(norm_list)

    ab_map_large = u_load(['Z:\asymmetry_project\data\mammograms\2004_screening\contralateral\abnormals\results\roi_' zerostr(fold, 2) '\' ab_list(ii).name]);
    ab_map = imresize(ab_map_large, 0.5, 'bilinear'); clear ab_map_large;
    norm_map_large = u_load(['Z:\asymmetry_project\data\mammograms\2004_screening\contralateral\normals\results\roi_' zerostr(fold, 2) '\' norm_list(ii).name]);
    norm_map = imresize(norm_map_large, 0.5, 'bilinear'); clear norm_map_large;
    
%     ab_image = u_load(['Z:\asymmetry_project\data\mammograms\2004_screening\contralateral\abnormals\' ab_list(ii).name(1:6)]);
%     norm_image = u_load(['Z:\asymmetry_project\data\mammograms\2004_screening\contralateral\normals\' norm_list(ii).name(1:6)]);

    min_prob = min(min(ab_map(:)), min(norm_map(:)));
    max_prob = max(max(ab_map(:)), max(norm_map(:)));
    
    if strcmpi(ab_list(ii).name(4), 'R')
        ab_map = fliplr(ab_map);
    else
        norm_map = fliplr(norm_map);
    end
    
    figure; 
%     subplot(2,2,1); imagesc(ab_image); axis image;
%     subplot(2,2,2); imagesc(norm_image); axis image;
    subplot(1,2,1); imagesc(1-ab_map); axis image; caxis([0 1]); colormap(jet(256));
    subplot(1,2,2); imagesc(1-norm_map); axis image; caxis([0 1]);
    clear ab_map norm_map;
end
%%
ab_sum = 0;
ab_N = 0;
norm_sum = 0;
norm_N = 0;
ab_means = [];
norm_means = [];

for fold = 1:3
    
    ab_list = dir(['Z:\asymmetry_project\data\mammograms\2004_screening\contralateral\abnormals\results\roi_' zerostr(fold, 2) '\*.mat']);
    norm_list = dir(['Z:\asymmetry_project\data\mammograms\2004_screening\contralateral\normals\results\roi_' zerostr(fold, 2) '\*.mat']);
    
    for ii = 1:length(norm_list)
        
        ab_map = u_load(['Z:\asymmetry_project\data\mammograms\2004_screening\contralateral\abnormals\results\roi_' zerostr(fold, 2) '\' ab_list(ii).name]);
        norm_map = u_load(['Z:\asymmetry_project\data\mammograms\2004_screening\contralateral\normals\results\roi_' zerostr(fold, 2) '\' norm_list(ii).name]);
        
        if strcmpi(ab_list(ii).name(4), 'R')
            ab_map = fliplr(ab_map);
        else
            norm_map = fliplr(norm_map);
        end
        ab_mask = u_load(['Z:\asymmetry_project\data\masks\2004_screening\contralateral\abnormals\' ab_list(ii).name(1:6) '_mask.mat']);
        norm_mask = u_load(['Z:\asymmetry_project\data\masks\2004_screening\contralateral\normals\' norm_list(ii).name(1:6) '_mask.mat']);       
        
        ab_im_sum = sum(1 - ab_map(ab_mask));
        ab_im_N = sum(ab_mask(:));
        norm_im_sum = sum(1 - norm_map(norm_mask));
        norm_im_N = sum(norm_mask(:));
        
        %display(['Image ' num2str(ii) ': ' num2str([ab_im_sum/ab_im_N norm_im_sum/norm_im_N], 3) ]); 
        
        ab_sum = ab_sum + ab_im_sum;
        ab_N = ab_N + ab_im_N;
        norm_sum = norm_sum + norm_im_sum;
        norm_N = norm_N + norm_im_N;
        
        ab_means(end+1) = ab_im_sum/ab_im_N; %#ok
        norm_means(end+1) = norm_im_sum/norm_im_N; %#ok
        
    end
end
display(ab_sum / ab_N);
display(norm_sum / norm_N);
[h p] = ttest(ab_means, norm_means) %#ok
sum(ab_means > norm_means) %#ok
%%
for fold = 3%1:3
    
    ab_list = dir(['Z:\asymmetry_project\data\mammograms\2004_screening\contralateral\abnormals\results\roi_' zerostr(fold, 2) '\*.mat']);
    norm_list = dir(['Z:\asymmetry_project\data\mammograms\2004_screening\contralateral\normals\results\roi_' zerostr(fold, 2) '\*.mat']);
    
    for ii = 8:length(norm_list)
        
        ab_map = u_load(['Z:\asymmetry_project\data\mammograms\2004_screening\contralateral\abnormals\results\roi_' zerostr(fold, 2) '\' ab_list(ii).name]);
        norm_map = u_load(['Z:\asymmetry_project\data\mammograms\2004_screening\contralateral\normals\results\roi_' zerostr(fold, 2) '\' norm_list(ii).name]);
        
        ab_roi1 = u_load(['Z:\asymmetry_project\data\mammograms\2004_screening\contralateral_roi\abnormals\results\roi_' zerostr(fold, 2) '\' ab_list(ii).name]);
        norm_roi1 = u_load(['Z:\asymmetry_project\data\mammograms\2004_screening\contralateral_roi\normals\results\roi_' zerostr(fold, 2) '\' norm_list(ii).name]);
        
        if strcmpi(ab_list(ii).name(4), 'R')
            ab_map = fliplr(ab_map);
            ab_roi1 = fliplr(ab_roi1);
        else
            norm_map = fliplr(norm_map);
            norm_roi1 = fliplr(norm_roi1);
        end
        load(['C:\isbe\asymmetry_project\data\mammograms\2004_screening\contralateral_abnormal_roi\' ab_list(ii).name(1:6)]);
        
        contralateral_pair.abnormal_pos = round(contralateral_pair.abnormal_pos/2);
        contralateral_pair.normal_pos = round(contralateral_pair.normal_pos/2);
        
        ra1 = contralateral_pair.abnormal_pos(1,2);
        ra2 = contralateral_pair.abnormal_pos(2,2);
        ca1 = contralateral_pair.abnormal_pos(1,1);
        ca2 = contralateral_pair.abnormal_pos(2,1);
        
        rn1 = contralateral_pair.normal_pos(1,2);
        rn2 = contralateral_pair.normal_pos(2,2);
        cn1 = contralateral_pair.normal_pos(1,1);
        cn2 = contralateral_pair.normal_pos(2,1);
        
        clear contralateral_pair;
        
        ab_roi2 = ab_map(ra1:ra2, ca1:ca2); clear ab_map;
        norm_roi2 = norm_map(rn1:rn2, cn1:cn2); clear norm_map;
        
        figure; 
        subplot(2,2,1); imagesc(ab_roi1); axis image; colormap(jet(256));
        subplot(2,2,2); imagesc(norm_roi1); axis image;
        subplot(2,2,3); imagesc(ab_roi2); axis image;
        subplot(2,2,4); imagesc(norm_roi2); axis image;
        
        %display([mean..
        
    end
end
%%
rel_list = dir('Z:\asymmetry_project\data\relevance_maps\2004_screening\abnormals\*.mat');
for ii = 1:length(rel_list)
    load(['Z:\asymmetry_project\data\relevance_maps\2004_screening\abnormals\' rel_list(ii).name]);
    probability_image = 1 - probability_image;
    save(['Z:\asymmetry_project\data\relevance_maps\2004_screening\abnormals\' rel_list(ii).name], 'probability_image');
    save(['C:\isbe\asymmetry_project\data\relevance_maps\2004_screening\abnormals\' rel_list(ii).name], 'probability_image');
    clear probability_image;
end
%%
for fold = 4:10
    
    ab_list = dir(['Z:\asymmetry_project\data\mammograms\2004_screening\contralateral\abnormals\results\roi_' zerostr(fold, 2) '\*.mat']);
    
    for ii = 1:length(ab_list)
        load(['Z:\asymmetry_project\data\mammograms\2004_screening\contralateral\abnormals\results\roi_' zerostr(fold, 2) '\' ab_list(ii).name]);
        probability_image = 1 - probability_image;
        
        save(['Z:\asymmetry_project\data\relevance_maps\2004_screening\abnormals\' ab_list(ii).name], 'probability_image');
        %save(['C:\isbe\asymmetry_project\data\relevance_maps\2004_screening\abnormals\' ab_list(ii).name], 'probability_image');
        clear probability_image;
    end
    
    norm_list = dir(['Z:\asymmetry_project\data\mammograms\2004_screening\contralateral\normals\results\roi_' zerostr(fold, 2) '\*.mat']);
    for ii = 1:length(norm_list)
        load(['Z:\asymmetry_project\data\mammograms\2004_screening\contralateral\normals\results\roi_' zerostr(fold, 2) '\' norm_list(ii).name]);
        probability_image = 1 - probability_image;
        
        save(['Z:\asymmetry_project\data\relevance_maps\2004_screening\abnormals\' norm_list(ii).name], 'probability_image');
        %save(['C:\isbe\asymmetry_project\data\relevance_maps\2004_screening\abnormals\' norm_list(ii).name], 'probability_image');
        clear probability_image;
    end
end
%%
r_list = dir('Z:\asymmetry_project\data\relevance_maps\2004_screening\abnormals\*R*.mat');
for ii = 1:length(r_list)
    load(['Z:\asymmetry_project\data\relevance_maps\2004_screening\abnormals\' r_list(ii).name]);
    probability_image = fliplr(probability_image);
    save(['Z:\asymmetry_project\data\relevance_maps\2004_screening\abnormals\' r_list(ii).name], 'probability_image');
end
copyfile('Z:\asymmetry_project\data\relevance_maps\2004_screening\abnormals\', 'C:\isbe\asymmetry_project\data\relevance_maps\2004_screening\abnormals\');
%%


%%
g_list = [...
    dir('Z:\asymmetry_project\data\radial_maps\2004_screening\abnormals\*064.mat');...
    dir('Z:\asymmetry_project\data\radial_maps\2004_screening\abnormals\*128.mat');...
    dir('Z:\asymmetry_project\data\radial_maps\2004_screening\abnormals\*256.mat')
    dir('Z:\asymmetry_project\data\radial_maps\2004_screening\abnormals\*inf.mat')];

for ii = 1:length(g_list)
    movefile(...
        ['Z:\asymmetry_project\data\radial_maps\2004_screening\abnormals\' g_list(ii).name],...
        ['Z:\asymmetry_project\data\weighted_radial_maps\2004_screening\abnormals\' g_list(ii).name]);
    copyfile(...
        ['Z:\asymmetry_project\data\weighted_radial_maps\2004_screening\abnormals\' g_list(ii).name],...
        ['G:\asymmetry_project\data\weighted_radial_maps\2004_screening\abnormals\' g_list(ii).name]);
end
%%
g_list = [...
    dir('Z:\asymmetry_project\data\radial_maps\2004_screening\normals\*128.mat');...
    dir('Z:\asymmetry_project\data\radial_maps\2004_screening\normals\*256.mat');...
    dir('Z:\asymmetry_project\data\radial_maps\2004_screening\normals\*inf.mat')];
    
for ii = 1:length(g_list)
    copyfile(...
        ['Z:\asymmetry_project\data\radial_maps\2004_screening\normals\' g_list(ii).name],...
        ['C:\isbe\asymmetry_project\data\radial_maps\2004_screening\normals\' g_list(ii).name]);
    copyfile(...
        ['Z:\asymmetry_project\data\radial_maps\2004_screening\normals\' g_list(ii).name],...
        ['G:\asymmetry_project\data\radial_maps\2004_screening\normals\' g_list(ii).name]);
end
%%
r_list = dir('Z:\asymmetry_project\data\relevance_maps\2004_screening\abnormals\*R*.mat');

for ii = 1:20;
    %for right breast
    r_map = u_load(['Z:\asymmetry_project\data\relevance_maps\2004_screening\abnormals\' r_list(ii).name]);
    l_map = u_load(['Z:\asymmetry_project\data\line_maps\2004_screening\abnormals\' r_list(ii).name(1:6) '_data.mat']);
    
    r_map = imresize(r_map, 0.5, 'bilinear');
    l_map = imresize(l_map, 0.5, 'bilinear');
    
    r_map = fliplr(r_map);
    write_im_from_colormap(l_map, ['G:\asymmetry_project\data\misc\figures\' r_list(ii).name(1:6) '_line_map.bmp'], gray(256), [0 1]);
    write_im_from_colormap(r_map .* l_map, ['G:\asymmetry_project\data\misc\figures\' r_list(ii).name(1:6) '_weighted_map.bmp'], gray(256), [0 1]);
    clear r_map l_map;
    
    l_name = r_list(ii).name;
    l_name(4) = 'L';
    
    r_map = u_load(['Z:\asymmetry_project\data\relevance_maps\2004_screening\abnormals\' l_name]);
    l_map = u_load(['Z:\asymmetry_project\data\line_maps\2004_screening\abnormals\' l_name(1:6) '_data.mat']);

    r_map = imresize(r_map, 0.5, 'bilinear');
    l_map = imresize(l_map, 0.5, 'bilinear');
    
    write_im_from_colormap(l_map, ['G:\asymmetry_project\data\misc\figures\' l_name(1:6) '_line_map.bmp'], gray(256), [0 1]);
    write_im_from_colormap(r_map .* l_map, ['G:\asymmetry_project\data\misc\figures\' l_name(1:6) '_weighted_map.bmp'], gray(256), [0 1]);
    
    clear r_map l_map;
end
%%
compute_rad_maps_batch('2004_screening\abnormals', 'relevance_dir', 'relevance_maps', 'view', 'LCC', 'num_jobs', 71, 'task_id', 24, 'radial_dir', 'weighted_radial_maps', 'distance_ranges', [64 128 256 inf]);
compute_rad_maps_batch('2004_screening\abnormals', 'relevance_dir', 'relevance_maps', 'view', 'LML', 'num_jobs', 75, 'task_id', 24, 'radial_dir', 'weighted_radial_maps', 'distance_ranges', [64 128 256 inf]);
%%