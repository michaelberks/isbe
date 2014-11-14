cd C:\isbe\matlab_code\mab\density\stepwedge\private\commented' code'\
%%
reduction_factor = 0.25;
patch = imresize(u_load('patch3.mat'), reduction_factor);
im_dir = 'C:\isbe\density\project2013\DigBatch5_anon\mammograms\';
%%
im_list = dir([im_dir '*RML*2430*.tif']);
left_breast = 0;
for i_im = 1:5
    mammo = rot90(imread([im_dir im_list(i_im).name]),2*left_breast);
    [markers, region] = find_markers_mb(mammo, patch, reduction_factor, left_breast, 1, 0);
end
%%
im_list = dir([im_dir '*LML*2430*.tif']);
left_breast = 1;
for i_im = 1:5
    mammo = rot90(imread([im_dir im_list(i_im).name]),2*left_breast);
    [markers, region] = find_markers_mb(mammo, patch, reduction_factor, left_breast, 1, 0);
end
%%
im_list = dir([im_dir '*RML*1824*.tif']);
left_breast = 0;
for i_im = 1:7
    mammo = rot90(imread([im_dir im_list(i_im).name]),1-2*left_breast);
    [markers, region] = find_markers_mb(mammo, patch, reduction_factor, left_breast, 1, 0);
end
%%
im_list = dir([im_dir '*LML*1824*.tif']);
left_breast = 1;
for i_im = 1:7
    mammo = rot90(imread([im_dir im_list(i_im).name]),1-2*left_breast);
    [markers, region] = find_markers_mb(mammo, patch, reduction_factor, left_breast, 1, 0);
end
%%
cd C:\isbe\matlab_code\mab\density\stepwedge\private\
%%
im_list = dir([im_dir '*RML*2430*.tif']);
left_breast = 0;
max_pairs = 4;
large_mammo = 1;

for i_im = 1:5
    mammo = rot90(imread([im_dir im_list(i_im).name]),2*left_breast);
    [x_pts, y_pts, selected_markers, marker_ui_xy] = ...
        marker_auto_detect(mammo, max_pairs, patch, large_mammo, left_breast, reduction_factor, 0);
end
%%
im_list = dir([im_dir '*LML*2430*.tif']);
left_breast = 1;
max_pairs = 4;
large_mammo = 1;

for i_im = 1:5
    mammo = rot90(imread([im_dir im_list(i_im).name]),2*left_breast);
    [x_pts, y_pts, selected_markers, marker_ui_xy] = ...
        marker_auto_detect(mammo, max_pairs, patch, large_mammo, left_breast, reduction_factor, 0);
end
%%
im_list = dir([im_dir '*RML*1824*.tif']);
left_breast = 0;
max_pairs = 3;
large_mammo = 0;

for i_im = 1:5
    mammo = rot90(imread([im_dir im_list(i_im).name]),1-2*left_breast);
    [x_pts, y_pts, selected_markers, marker_ui_xy] = ...
        marker_auto_detect(mammo, max_pairs, patch, large_mammo, left_breast, reduction_factor, 0);
end
%%
im_list = dir([im_dir '*LML*1824*.tif']);
left_breast = 1;
max_pairs = 3;
large_mammo = 0;

for i_im = 1
    mammo = rot90(imread([im_dir im_list(i_im).name]),1-2*left_breast);
    [x_pts, y_pts, selected_markers, marker_ui_xy] = ...
        marker_auto_detect(mammo, max_pairs, patch, large_mammo, left_breast, reduction_factor, 0);
end