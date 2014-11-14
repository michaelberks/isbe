mammo_dir = 'M:\asymmetry_project\data\mammograms\2004_screening\abnormals\mat\';
rf_dir = 'M:\asymmetry_project\data\contralateral_rfs\2004_screening\contralateral_images_binary\abnormals\';
maps_dir = 'M:\asymmetry_project\data\contralateral_rfs\2004_screening\abnormals\';
r_list = dir([mammo_dir '*RCC*.mat']);
l_list = dir([mammo_dir '*LCC*.mat']);
rf_list = dir([rf_dir 'random_forest*']);
%%
pack;
for ii = 62:71
    rf = u_load([rf_dir rf_list(ii).name]);
    class_map_r = rf.image1_votes1(2:2:end,2:2:end) ./ rf.image1_total_votes(2:2:end,2:2:end);
    class_map_r(~rf.image1_votes1(2:2:end,2:2:end)) = 0;
    
    class_map_l = rf.image2_votes1(2:2:end,2:2:end) ./ rf.image2_total_votes(2:2:end,2:2:end);
    class_map_l(~rf.image2_votes1(2:2:end,2:2:end)) = 0;
    class_map_l = fliplr(class_map_l);
    clear rf;
    
    figure;
    subplot(1,2,1); imagesc(class_map_r); axis image; colormap(jet(256));
    subplot(1,2,2); imagesc(class_map_l); axis image; colormap(jet(256));
    
    save([maps_dir r_list(ii).name(1:end-4) '_rf_map.mat'], 'class_map_r');
    save([maps_dir l_list(ii).name(1:end-4) '_rf_map.mat'], 'class_map_l');
    clear class_map_*;
end
    
