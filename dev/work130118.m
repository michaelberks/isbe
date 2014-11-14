data_dir = 'C:\isbe\asymmetry_project\data\retinograms\STARE\training\';
mkdir([data_dir 'vessel_masks\png\']);
for i_im = 1:20
    load([data_dir 'vessel_masks\' zerostr(i_im, 2) '_training_v_mask.mat']);
    imwrite(~v_mask, [data_dir 'vessel_masks\png\' zerostr(i_im, 2) '_training_v_mask.png']);
end
%%
data_dir = 'C:\isbe\asymmetry_project\data\fibre\training\';
f_list = dir([data_dir 'fibre_masks\*.mat']);
b_list = dir([data_dir 'fov_masks\*.mat']);

f_sum = 0;
b_sum = 0;

for i_im = 1:200
    f_mask = u_load([data_dir 'fibre_masks\' f_list(i_im).name]);
    b_mask = u_load([data_dir 'fov_masks\' b_list(i_im).name]);
    
    f_sum = f_sum + sum(f_mask(:) & b_mask(:));
    b_sum = b_sum + sum(~f_mask(:) & b_mask(:));
end

f_sum
b_sum