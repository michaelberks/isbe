plot_dir = 'K:\isbe\conferences_and_symposia\IPH_away_day\figures\';
%%
load('C:\isbe\asymmetry_project\data\retinograms\DRIVE\test\images\02_test.mat')
load('C:\isbe\asymmetry_project\data\retinograms\DRIVE\test\fov_masks\02_test_f_mask.mat');
load('C:\isbe\asymmetry_project\data\retinograms\DRIVE\test\vessel_masks\02_test_v_mask.mat');

ret_bg_r = ret(:,:,1);
ret_bg_r(vessel_mask | ~foveal_mask) = 255;
ret_bg_g = ret(:,:,2);
ret_bg_g(vessel_mask | ~foveal_mask) = 255;
ret_bg_b = ret(:,:,3);
ret_bg_b(vessel_mask | ~foveal_mask) = 255;

ret_fg_r = ret(:,:,1);
ret_fg_r(~vessel_mask) = 255;
ret_fg_g = ret(:,:,2);
ret_fg_g(~vessel_mask) = 255;
ret_fg_b = ret(:,:,3);
ret_fg_b(~vessel_mask) = 255;

ret_bg = cat(3, ret_bg_r, ret_bg_g, ret_bg_b);
ret_fg = cat(3, ret_fg_r, ret_fg_g, ret_fg_b);

figure; imgray(ret_bg);
figure; imgray(ret_fg);

imwrite(ret_fg, [plot_dir 'split_ret_fg.png']);
imwrite(ret_bg, [plot_dir 'split_ret_bg.png']);

bg_trans_map = double(~vessel_mask & foveal_mask);
fg_trans_map = double(vessel_mask);

imwrite(ret, [plot_dir 'split_ret_fg2.png'], 'alpha', fg_trans_map);
imwrite(ret, [plot_dir 'split_ret_bg2.png'], 'alpha', bg_trans_map);
%%
ret_rgb = rgb2gray(ret);
[bg_r bg_c] = find(bg_trans_map);
[fg_r fg_c] = find(fg_trans_map);
num_bg = length(bg_r);
num_fg = length(fg_r);
for ii = 1:20
    
    bg_i = ceil(rand*num_bg);
    fg_i = ceil(rand*num_fg);
    
    bg_patch = double(sample_window(ret_rgb, 33, bg_r(bg_i), bg_c(bg_i)));
    fg_patch = double(sample_window(ret_rgb, 33, fg_r(fg_i), fg_c(fg_i)));
    
    write_im_from_colormap(bg_patch, ...
        [plot_dir 'patches_bg ' zerostr(ii,2) '.png']);
    write_im_from_colormap(fg_patch, ...
        [plot_dir 'patches_fg ' zerostr(ii,2) '.png']);
end
%%
fg_ir = 449;
fg_ic = 360;
fg_patch = double(sample_window(ret_rgb, 33, fg_ir, fg_ic));
fg_patch_i = -interp2(double(fg_patch),1,'bicubic');

make_3d_patch_video(fg_patch_i, [plot_dir 'fg_patch_3d_vid.gif'], 'colormap', flipud(gray(256)));
make_improfile_video(fg_patch_i, [plot_dir 'fg_patch_ip_vid.gif'], 32, 30, 'colormap', flipud(gray(256)));
make_improfile_plots(fg_patch_i, [plot_dir 'fg_patch_ip_plots'], 32, 30, 'colormap', flipud(gray(256)));
%%
bg_ir = 439;
bg_ic = 350;
bg_patch = double(sample_window(ret_rgb, 33, bg_ir, bg_ic));
bg_patch_i = -interp2(double(bg_patch),1,'bicubic');
make_improfile_plots(bg_patch_i, [plot_dir 'bg_patch_ip_plots'], 32, 30, 'colormap', flipud(gray(256)));
%%
od_ir = 305;
od_ic = 425;
od_patch = double(sample_window(ret_rgb, 33, od_ir, od_ic));
od_patch_i = -interp2(double(od_patch),1,'bicubic');

make_3d_patch_video(od_patch_i, [plot_dir 'od_patch_3d_vid.gif'], 'colormap', flipud(gray(256)));
make_improfile_video(od_patch_i, [plot_dir 'od_patch_ip_vid.gif'], 33, 33, 'colormap', flipud(gray(256)));
make_improfile_plots(od_patch_i, [plot_dir 'od_patch_ip_plots'], 33, 33, 'colormap', flipud(gray(256)));
%%
load(['C:\isbe\asymmetry_project\data\retinograms\DRIVE\test\images\' zerostr(14,2) '_test.mat']);
ret_rgb = rgb2gray(ret);
ha_ir = 305;
ha_ic = 248;
ha_patch = double(sample_window(ret_rgb, 33, ha_ir, ha_ic));
ha_patch_i = -interp2(double(ha_patch),1,'bicubic');
make_improfile_plots(ha_patch_i, [plot_dir 'ha_patch_ip_plots'], 33, 33, 'colormap', flipud(gray(256)));
%%
ve_ir = 317;
ve_ic = 151;
ve_patch = double(sample_window(ret_rgb, 33, ve_ir, ve_ic));
ve_patch_i = -interp2(double(ve_patch),1,'bicubic');
make_improfile_plots(ve_patch_i, [plot_dir 've_patch_ip_plots'], 33, 33, 'colormap', flipud(gray(256)));
%%
[g dg ddg] = gaussian_filters_1d(8);
g = [zeros(1,40) g zeros(1,40)];
hg = g;
hg(ceil(end/2):end) = hg(ceil(end/2));

figure; plot(g, 'r', 'linewidth', 4); axis off; exportfig([plot_dir 'profile_line.png']);
figure; plot(hg, 'r', 'linewidth', 4); axis off; exportfig([plot_dir 'profile_edge.png']);
figure; plot(dg, 'b', 'linewidth', 4); axis off; exportfig([plot_dir 'filter_odd_1d.png']);
figure; plot(ddg, 'b',  'linewidth', 4); axis off; exportfig([plot_dir 'filter_even_1d.png']);
%
figure; hold on;
plot(conv(g, dg, 'valid'), 'm', 'linewidth', 4);  
plot(get(gca, 'xlim'), [0 0], 'k--', 'linewidth', 2);
plot([41 41], get(gca, 'ylim'), 'k--', 'linewidth', 2);
axis off; exportfig([plot_dir 'response_line_odd_1d.png']);

figure; hold on;
plot(conv(g, ddg, 'valid'), 'm', 'linewidth', 4);
plot(get(gca, 'xlim'), [0 0], 'k--', 'linewidth', 2);
plot([41 41], get(gca, 'ylim'), 'k--', 'linewidth', 2);
axis off; exportfig([plot_dir 'response_line_even_1d.png']);
 
figure; hold on;
plot(conv(hg, dg, 'valid'), 'm', 'linewidth', 4); axis off; 
plot(get(gca, 'xlim'), [0 0], 'k--', 'linewidth', 2);
plot([31.5 31.5], get(gca, 'ylim'), 'k--', 'linewidth', 2);
exportfig([plot_dir 'response_edge_odd_1d.png']);

figure; hold on;
plot(conv(hg, ddg, 'valid'), 'm', 'linewidth', 4);
plot(get(gca, 'xlim'), [0 0], 'k--', 'linewidth', 2);
plot([31.5 31.5], get(gca, 'ylim'), 'k--', 'linewidth', 2);
axis off; exportfig([plot_dir 'response_edge_even_1d.png']);
%%
%---------------------------------------------
giant_patch = imread('C:\isbe\nailfold\data\rsa_study\apexes\giant\apex0017.png');
load('C:\isbe\nailfold\data\rsa_study\apexes\normal\apex2131.mat');
apex_xy(:,2) = apex_xy(:,2) + apex_properties.sr;
apex_xy(:,1) = apex_xy(:,1) + apex_properties.sc;
apex_xyc = mean(apex_xy);
normal_patch = sample_window(nailfold, 401, round(apex_xyc(2)), round(apex_xyc(1)));

write_im_from_colormap(giant_patch, [plot_dir 'giant_patch.png']);
write_im_from_colormap(normal_patch, [plot_dir 'normal_patch.png']);

make_improfile_plots(-double(normal_patch), [plot_dir 'nailfold_normal_ip_plots'],...
    194, 231, 'colormap', flipud(gray(256)), 'ylims', [-180 -150]);

make_improfile_plots(-double(giant_patch), [plot_dir 'nailfold_giant_ip_plots'],...
    96, 264, 'colormap', flipud(gray(256)), 'ylims', [-150 -90], 'profile_radius', 20);

make_improfile_plots(-double(giant_patch), [plot_dir 'nailfold_giant_ip_plots2'],...
    96, 264, 'colormap', flipud(gray(256)), 'ylims', [-150 -90], 'profile_radius', 60);
%%
%--------------------------------------------------------------------------
image_dir = 'C:\isbe\nailfold\data\rsa_study\training\images\';
pred_dir = 'C:\isbe\nailfold\data\rsa_study\training\predictions\detection\rf_classification\182321\';
ori_dir = 'C:\isbe\nailfold\data\rsa_study\training\predictions\orientation\rf_regression\182263\';
width_dir = 'C:\isbe\nailfold\data\rsa_study\training\predictions\width\rf_regression\182367\';

i_list = dir([image_dir '*.mat']);
p_list = dir([pred_dir '*.mat']);
o_list = dir([ori_dir '*.mat']);
w_list = dir([width_dir '*.mat']);

i_ve = 23;
vessel_im = u_load([image_dir i_list(i_ve).name]);
vessel_pred = u_load([pred_dir p_list(i_ve).name]);
vessel_width = u_load([width_dir w_list(i_ve).name]);
vessel_ori = u_load([ori_dir o_list(i_ve).name]);

ori_im = complex2rgb(vessel_pred.*vessel_ori);
figure; 
subplot(1,2,1); imgray(vessel_im);
subplot(1,2,2); imgray(two_map_colour(vessel_pred, vessel_width, [], [10 20]));

figure; 
subplot(1,2,1); imgray(vessel_im);
subplot(1,2,2); imgray();

write_im_from_colormap(vessel_im(66:end,:), [plot_dir 'nailfold_vessels.png']);
imwrite(ori_im(66:end,:,:), [plot_dir 'nailfold_vessels_ori.png']);
%%
point_src = zeros(128);
point_src(64,64) = 1;
[h2d_responses] = compute_hilbert_2nd_derivatives_sep(point_src, 16);
h2d_f0 = steer_hilbert_2nd_derivatives(h2d_responses, pi/2);
h2d_f1 = steer_hilbert_2nd_derivatives(h2d_responses, pi/6);
h2d_f2 = steer_hilbert_2nd_derivatives(h2d_responses, -pi/6);

figure;
subplot(1,3,1); imgray(h2d_f0);
subplot(1,3,2); imgray(h2d_f1);
subplot(1,3,3); imgray(h2d_f2);
%%
write_im_from_colormap(h2d_f0, [plot_dir 'h2d_0.png']);
write_im_from_colormap(h2d_f1, [plot_dir 'h2d_1.png']);
write_im_from_colormap(h2d_f2, [plot_dir 'h2d_2.png']);

[h2d_responses] = compute_hilbert_2nd_derivatives_sep(point_src, 8);
h2d_f0 = steer_hilbert_2nd_derivatives(h2d_responses, pi/2);
write_im_from_colormap(h2d_f0, [plot_dir 'h2d_0_1.png']);

[h2d_responses] = compute_hilbert_2nd_derivatives_sep(point_src, 4);
h2d_f0 = steer_hilbert_2nd_derivatives(h2d_responses, pi/2);
write_im_from_colormap(h2d_f0, [plot_dir 'h2d_0_2.png']);
    
    