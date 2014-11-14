data_dir = 'C:\isbe\asymmetry_project\data\retinograms\STARE\training\';
im_list = dir([data_dir 'images\*.ppm']);
gt_list = dir([data_dir 'vessel_masks\ah\*.ppm']);

for i_im = 1:20
    ret = imread([data_dir 'images\' im_list(i_im).name]);
    v_mask = imread([data_dir 'vessel_masks\ah\' gt_list(i_im).name]);
    v_mask = v_mask > 0;    
    f_mask = ret(:,:,1) > 40;
      
    figure;
    subplot(1,2,1); imgray(ret);
    subplot(1,2,2); imgray(v_mask);
    
    save([data_dir 'images\' zerostr(i_im, 2) '_training.mat'], 'ret');
    save([data_dir 'vessel_masks\' zerostr(i_im, 2) '_training_v_mask.mat'], 'v_mask');
    save([data_dir 'fov_masks\' zerostr(i_im, 2) '_training_f_mask.mat'], 'f_mask');
end
%%
data_dir = 'C:\isbe\asymmetry_project\data\retinograms\STARE\training\';
gt_list = dir([data_dir 'vessel_masks\vk\*.ppm']);
for i_im = 1:20
    mask2 = imread([data_dir 'vessel_masks\vk\' gt_list(i_im).name]);
    mask2 = mask2 > 0;    
    mask = u_load([data_dir 'vessel_masks\' zerostr(i_im, 2) '_training_v_mask.mat']);
    fov_mask = u_load([data_dir 'fov_masks\' zerostr(i_im, 2) '_training_f_mask.mat']);
    
    t_pos(ii) = sum(mask2(mask & fov_mask));
    f_pos(ii) = sum(mask2(~mask & fov_mask));
    pos(ii) = sum(mask(:) & fov_mask(:));
    neg(ii) = sum(~mask(:) & fov_mask(:));
    
    temp = mask;
    mask = mask2;
    mask2 = temp;
    
    t_pos2(ii) = sum(mask2(mask & fov_mask));
    f_pos2(ii) = sum(mask2(~mask & fov_mask));
    pos2(ii) = sum(mask(:) & fov_mask(:));
    neg2(ii) = sum(~mask(:) & fov_mask(:));
    
    figure; 
    subplot(2,2,1); imgray(mask);
    subplot(2,2,2); imgray(mask2);
    subplot(2,2,3); imgray(mask & ~mask2);
    subplot(2,2,4); imgray(~mask & mask2);
end

display(sum(t_pos) / sum(pos))
display(sum(f_pos) / sum(neg))
%%
data_dir = 'C:\isbe\asymmetry_project\data\retinograms\STARE\training\';
mkdir([data_dir 'fov_masks_eroded']);
for i_im = 1:20
    f_mask = u_load([retroot 'fov_masks\' zerostr(i_im,2) '_training_f_mask.mat']);
    f_mask = imerode(f_mask, strel('disk', 10));
    save([data_dir 'fov_masks_eroded\' zerostr(i_im, 2) '_training_f_mask.mat'], 'f_mask');
end
%%
