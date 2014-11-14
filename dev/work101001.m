load('C:\isbe\asymmetry_project\data\synthetic_lines\lines512\image002.mat');
line_map = u_load('C:\isbe\asymmetry_project\data\synthetic_lines\lines512\results\191906\probability_image002.mat');
ori_map = u_load('C:\isbe\asymmetry_project\data\synthetic_lines\lines512\results\191934\probability_image002.mat');



[w_ori] = display_weighted_orientation(test_image, ori_map, line_map, 4, gray(256));
display_orientation(test_image, ori_map, 8, label > 0, gray(256));
for sigma = 2:10
    w_ori_smooth = imfilter(w_ori, fspecial('gaussian', 5*sigma, sigma));
    w_ori_map = 180*angle(w_ori_smooth)/pi;
    display_orientation(test_image, w_ori_map, 8, label > 0, gray(256));
end