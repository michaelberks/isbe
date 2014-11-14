%BII image comp script
vessel_patch = double(rgb2gray(imread('C:\isbe\nailfold\images\BII_comp\10598c_small.png')));

rf_ori = u_load('C:\isbe\nailfold\models\vessel\orientation\rf_regression\259076\predictor.mat');
rf_ori.tree_root = 'C:\isbe\nailfold\models\vessel\orientation\rf_regression/';
ori_args = u_load('C:\isbe\nailfold\models\vessel\orientation\rf_regression\259076\job_args.mat');

rf_det = u_load('C:\isbe\nailfold\models\vessel\detection\rf_classification\257273\predictor.mat');
rf_det.tree_root = 'C:\isbe\nailfold\models\vessel\detection\rf_classification/';
det_args = u_load('C:\isbe\nailfold\models\vessel\detection\rf_classification\257273\job_args.mat');

rf_wid = u_load('C:\isbe\nailfold\models\vessel\width\rf_regression\257847\predictor.mat');
rf_wid.tree_root = 'C:\isbe\nailfold\models\vessel\width\rf_regression/';
wid_args = u_load('C:\isbe\nailfold\models\vessel\width\rf_regression\257847\job_args.mat');

[patch_ori] = predict_image(...
    'image_in', vessel_patch,...
    'decomposition_args', ori_args.decomposition_args,...
    'predictor', rf_ori, ...
    'prediction_type', 'rf_regression',...
    'output_type', 'orientation',...
    'use_probs', 0,...
    'mask', [],...
    'tree_mask', [], ...
    'num_trees', [], ...
    'max_size', 128,...
    'incremental_results', 0);
%
[patch_det] = predict_image(...
    'image_in', vessel_patch,...
    'decomposition_args', det_args.decomposition_args,...
    'predictor', rf_det, ...
    'prediction_type', 'rf_classification',...
    'output_type', 'detection',...
    'use_probs', 0,...
    'mask', [],...
    'tree_mask', [], ...
    'num_trees', [], ...
    'max_size', 128,...
    'incremental_results', 0);
%
[patch_wid] = predict_image(...
    'image_in', vessel_patch,...
    'decomposition_args', wid_args.decomposition_args,...
    'predictor', rf_wid, ...
    'prediction_type', 'rf_regression',...
    'output_type', 'width',...
    'use_probs', 0,...
    'mask', [],...
    'tree_mask', [], ...
    'num_trees', [], ...
    'max_size', 128,...
    'incremental_results', 0);
%
figure; imgray(complex2rgb(patch_ori));
figure; imgray(patch_det);
figure; imgray(two_map_colour(patch_det, patch_wid, [], [15 50]));

save('C:\isbe\nailfold\images\BII_comp\10598c_small_ori.mat', 'patch_ori');
save('C:\isbe\nailfold\images\BII_comp\10598c_small_det.mat', 'patch_det');
save('C:\isbe\nailfold\images\BII_comp\10598c_small_wid.mat', 'patch_wid');
%%
load([nailfoldroot 'data\apex_detection\apex_templates_aligned.mat'], 'apex_template*');
load([nailfoldroot 'data\apex_detection\aligned\thresh.mat'], 'thresh2');
load([nailfoldroot 'data\apex_detection\mean_shape.mat'], 'mean_shape');

pilot_dir = 'C:\isbe\nailfold\data\pilot_study\';

aam_path = 'C:\isbe\nailfold\images\BII_comp\template_matching\aligned\';

detect_capillaries('C:\isbe\nailfold\images\BII_comp\10598c_small.png', ...
        'g2d_template', apex_template2_a,...
        'template_thresh', thresh2,...
        'max_num_candidates', 11,...
        'thetas', -15:3:15,...
        'scales', 0.8:0.1:1.5,...
        'mean_shape', mean_shape,...
        'aam_dir', aam_path,...
        'aam_exe', 'ncm_sandpit_mb',...
        'aam_path', 'C:\isbe\nailfold\playground\model_experiments\models\orig\2\vessel_apex_orig.smd',...
        'vessel_probs', [],...
        'prob_thresh', 0.9,...
        'do_template_matching', 1,...
        'do_aam', 1, ...
        'do_final_vessel', 1, ...
        'delete_candidate_patches', 0,...
        'nailfold_mask', []);
%%

figure; imgray(vessel_patch);
distal_apexes = zeros(0,2);
for i_ap = 1:11
    load([aam_path 'apex' zerostr(i_ap, 4) '_aam.mat'], 'apex_candidate');
    if isfield(apex_candidate, 'is_distal') && apex_candidate.is_distal
        apex_x = apex_candidate.fitted_vessel_xy(:,1) + apex_candidate.sc;
        apex_y = apex_candidate.fitted_vessel_xy(:,2) + apex_candidate.sr;
        plot(apex_x, apex_y, 'x');
        distal_apexes(end+1,:) = round([apex_x(16) apex_y(16)]);
    end

end
%%
[~,~, path_map] = shape_context_prob_track_mult(...
                    patch_ori.*patch_det, patch_det, patch_ori,...
                    distal_apexes(:,2), distal_apexes(:,1),...
                    'num_streams', 1e4, 'stopping_prob', 0.1, 'min_rdist', 8);
%
for i_ap = 1:11
    if isfield(apex_candidate, 'is_distal') && apex_candidate.is_distal
        path_map(distal_apexes(i_ap,2), distal_apexes(i_ap,1)) = 501;
    end

end
complex_path_map = (path_map > 500) .* exp(1i*angle(patch_ori));
path_map_rgb = complex2rgb(complex_path_map);

vessel_im = (vessel_patch - min(vessel_patch(:))) / (max(vessel_patch(:))-min(vessel_patch(:)));
vessel_im = cat(3, vessel_im, vessel_im, vessel_im);
%
alpha = 0.75;

vessel_path_rgb = alpha*vessel_im + (1-alpha)*path_map_rgb;



%%
alpha = 0.8;
vessel_thresh = bwselect(patch_det > 0.8, distal_apexes(:,1) , distal_apexes(:,2));
complex_det_map = vessel_thresh .* patch_det .* exp(1i*angle(patch_ori));
det_map_rgb = complex2rgb(complex_det_map);
vessel_det_rgb = alpha*vessel_im + (1-alpha)*det_map_rgb;
%%
[nrows, ncols] = size(vessel_patch);
alpha = repmat(1-linspace(0,1,ncols).^3, nrows, 1);
alpha = cat(3, alpha, alpha, alpha);
vessel_det_fade = alpha.*vessel_im + (1-alpha).*det_map_rgb;

%%
figure; imgray(vessel_patch);
figure; imgray(complex2rgb(patch_ori));
figure; imgray(patch_det);

figure; imgray(det_map_rgb);
figure; imgray(path_map_rgb);

figure; imgray(vessel_det_rgb);
figure; imgray(vessel_det_fade);
figure; imgray(vessel_path_rgb);
%%
write_im_from_colormap(vessel_patch, 'C:\isbe\nailfold\images\BII_comp\images\00_nailfold.png');
imwrite(complex2rgb(patch_ori), 'C:\isbe\nailfold\images\BII_comp\images\orientation_prediction.png');
imwrite(patch_det, 'C:\isbe\nailfold\images\BII_comp\images\vessel_prediction.png');

imwrite(det_map_rgb, 'C:\isbe\nailfold\images\BII_comp\images\vessel_masked_orientation.png');
imwrite(path_map_rgb, 'C:\isbe\nailfold\images\BII_comp\images\apex_stream_orientation.png');

imwrite(vessel_det_rgb, 'C:\isbe\nailfold\images\BII_comp\images\vessel_orientation_overlay.png');
imwrite(vessel_det_fade, 'C:\isbe\nailfold\images\BII_comp\images\vessel_orientation_fade.png');
imwrite(vessel_path_rgb, 'C:\isbe\nailfold\images\BII_comp\images\apex_stream_orientation_overlay.png');

%%
apex_radius = 5;
max_apex_guess = 500;
vessel_prob_smoothing_sigma = 2;
curvature_smoothing_sigma = 0.5;
strong_vessel_thresh = 0.25;
curv_max = 0.5;
smoothing_window = ones(1,25) / 25;

g = gaussian_filters_1d(vessel_prob_smoothing_sigma);
g = g / sum(g);

vessel_ori = u_load('C:\isbe\nailfold\images\BII_comp\10598c_small_ori.mat');
vessel_prob = u_load('C:\isbe\nailfold\images\BII_comp\10598c_small_det.mat');

vessel_prob = conv2(g', g, vessel_prob, 'same');
%vessel_width = conv2(g', g, vessel_width, 'same');
    
%Compute NMS centrelines
vessel_nms = mb_non_maximal_supp(vessel_prob, angle(vessel_ori)/2);
strong_vessels = vessel_nms > strong_vessel_thresh;
if any(strong_vessels(:))
    [rstrong cstrong] = find(strong_vessels);
    vessel_centre_mask = bwselect(vessel_nms > 0, cstrong, rstrong, 8);
else
    vessel_centre_mask = strong_vessels;
end
[vessel_centre_y vessel_centre_x] = find(vessel_centre_mask);

figure; imgray(vessel_centre_mask);

imwrite(vessel_centre_mask, 'C:\isbe\nailfold\images\BII_comp\images\vessel_centre.png');
%%
[nrows ncols] = size(vessel_centre_mask);
f1 = figure('windowstyle', 'normal', 'position', [100 100 ncols nrows]); 
a1 = axes('units', 'pixels', 'position', [0 0 ncols nrows]);
imgray(vessel_centre_mask);
%
vxy_pts = [450 159 586; 128 308 80];
x = repmat(-24:0.5:24, length(-24:0.5:24), 1);
y = x';
xy = [x(:) y(:)];

colors = 'rgb';
plot_num = 1;
for vxy = vxy_pts
    
    ori_c = angle(vessel_ori(vxy(2), vxy(1)))/2;
    width_c = patch_wid(vxy(2), vxy(1));
    
    rot = [cos(ori_c) -sin(ori_c); sin(ori_c) cos(ori_c)];
    scale = width_c / 20;

    %Transform points given scale and angle and translate to
    %candidate position
    xya = xy * rot * scale;
    xa = reshape(xya(:,1) + vxy(1), size(x));
    ya = reshape(xya(:,2) + vxy(2), size(y));
    
    plot(a1,...
        [xa(1,1) xa(1,end) xa(end,end) xa(end,1) xa(1,1)],...
        [ya(1,1) ya(1,end) ya(end,end) ya(end,1) ya(1,1)], colors(plot_num), 'linewidth', 3);
    plot_num = plot_num + 1;
    
    vessel_prob_patch = interp2(vessel_prob, xa, ya);
    
    figure; imgray(vessel_prob_patch);
    imwrite(vessel_prob_patch, ['C:\isbe\nailfold\images\BII_comp\images\vessel_centre_x'...
        zerostr(vxy(1),3) 'y' zerostr(vxy(2),3) '.png']);
end
figure(f1);
exportfig('C:\isbe\nailfold\images\BII_comp\images\vessel_centre_boxes.png');
%%
load('C:\isbe\nailfold\images\BII_comp\10598c_small_ori.mat', 'patch_ori');
load('C:\isbe\nailfold\images\BII_comp\10598c_small_det.mat', 'patch_det');
load('C:\isbe\nailfold\images\BII_comp\10598c_small_wid.mat', 'patch_wid');

apex_class_rf = u_load('C:\isbe\nailfold\models\apex/classification/frog/rf.mat');
apex_offset_x_rf = u_load('C:\isbe\nailfold\models\apex/offset_x/frog/rf.mat');
apex_offset_y_rf = u_load('C:\isbe\nailfold\models\apex/offset_y/frog/rf.mat');

g = gaussian_filters_1d(2);
g = g / sum(g);
%

hog_args.cell_sz = [8 8];
hog_args.num_ori_bins = 9;

[vessel_centre] = extract_vessel_centres(patch_det, patch_ori, patch_wid);

[apex_offset_map apex_class_pred apex_offset_x_pred apex_offset_y_pred] = ...
        predict_apex_offsets(...
            'apex_class_rf', apex_class_rf,...
            'apex_offset_x_rf', apex_offset_x_rf,...
            'apex_offset_y_rf', apex_offset_y_rf,...
            'vessel_feature_im', patch_det, ...
            'vessel_centre', vessel_centre, ...
            'smoothing_sigma', 2,...
            'num_cells', 8,...
            'hog_args', hog_args,...
            'apex_class_thresh', 0.5);
        
[nrows ncols] = size(nrows, ncols);
apex_class_map = zeros(nrows, ncols);
vessel_centre_idx = sub2ind(nrows, ncols, vessel_centre.y, vessel_centre.x);
apex_class_map(vessel_centre_idx) =  apex_class_pred;

write_im_from_colormap(apex_class_map, 'C:\isbe\nailfold\images\BII_comp\images\apex_class_map.png', gray(256));
write_im_from_colormap(apex_offset_map, 'C:\isbe\nailfold\images\BII_comp\images\apex_offset_map.png', hot(256));

%%
% [apex_offset_map apex_class_pred apex_offset_x_pred apex_offset_y_pred] = ...
%     predict_apex_offsets(apex_class_rf, apex_offset_x_rf, apex_offset_y_rf,...
%         vessel_feature_im, vessel_centre, hog_args, xy, hog_sz, 0.5, 1000, 20);

vessel_patch = double(imread('C:\isbe\nailfold\images\BII_comp\10598c_small.png'));
vessel_im = (vessel_patch - min(vessel_patch(:))) / (max(vessel_patch(:))-min(vessel_patch(:)));

apex_offset_im = ind2rgb(ceil(256*apex_offset_patch/max(apex_offset_patch(:))), hot(256));

alpha = 0.75;
vessel_offset_rgb = alpha*vessel_im + (1-alpha)*apex_offset_im;
imwrite(vessel_offset_rgb, 'C:\isbe\nailfold\images\BII_comp\images\vessel_offset_rgb.png');