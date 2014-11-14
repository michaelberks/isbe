rsa_dir = 'rsa_study/';
image_dir = [nailfoldroot 'data/' rsa_dir 'test/images/'];
fov_mask_dir = [nailfoldroot 'data/' rsa_dir 'test/fov_masks/'];
prob_dir = [nailfoldroot 'data/' rsa_dir 'test/predictions/detection/rf_classification/222836/'];
ori_dir = [nailfoldroot 'data/' rsa_dir 'test/predictions/orientation/rf_regression/222835/'];

test_list = dir([image_dir '*.mat']);
num_images = length(test_list);
%
do_dt = 0;
for i_im = 21:40%num_images;
    
    im_num = test_list(i_im).name(1:6);
    display(['Processing image ' num2str(i_im) ' of ' num2str(num_images)]);

    nailfold = u_load([image_dir im_num '.mat']);
    f_mask = u_load([fov_mask_dir im_num '_f_mask.mat']);
    f_mask = imerode(f_mask, strel('disk', 50));
    nailfold_prob = u_load([prob_dir im_num '_pred.mat']);
    nailfold_ori = u_load([ori_dir im_num '_pred.mat']);
    
    nailfold = imresize(nailfold, 0.5);
    f_mask = imresize(f_mask, 0.5);
    nailfold_prob = imresize(nailfold_prob, 0.5, 'bilinear');
    nailfold_ori = imresize(nailfold_ori, 0.5, 'bilinear');
    
    nailfold_prob(~f_mask) = 0;
    
    num_pts = sum(f_mask(:));
    
    if do_dt
        num_levels = 5;
        nailfold_dt = compute_dual_tree(nailfold, num_levels, 0);

        %
        nailfold_vec = zeros(6, 5);
        for i_levels = 1:5
            f_mask_lev = imresize(f_mask, 1/(2^i_levels));  
            nailfold_prob_lev = imresize(nailfold_prob, 1/(2^i_levels));  
            for i_bands = 1:6

                subband = nailfold_dt{i_levels}(:,:,i_bands);
                %nailfold_vec(i_bands, i_levels) = sum(abs(subband(f_mask_lev)));
                nailfold_vec(i_bands, i_levels) = sum(abs(subband(:) .* nailfold_prob_lev(:)));
            end
        end
    else
        vessel_curv = -(complex_gaussian_curvature(nailfold_ori ./ (abs(nailfold_ori)+1e-6), 0.5));

        nailfold_vec = compute_weighted_histogram(vessel_curv, nailfold_prob, linspace(-0.1, 0.1, 100));
    end
    
    nailfold_vec = nailfold_vec / num_pts;   
    
    figure; 
    subplot(2,1,1);
    imgray(nailfold); caxis([min(nailfold(f_mask)) max(nailfold(f_mask))]);
    
    subplot(2,1,2);
    bar(1:length(nailfold_vec(:)), nailfold_vec(:));
end
%%
load('C:\isbe\nailfold\models\vessel\detection\rf_classification\257273\01_sampled_data\image_lists.mat');
[image_list{1:length(image_lists),1}] = image_lists.image;
curr_tree = 1;
sampled_maps_dir = 'C:\isbe\nailfold\models\vessel\detection\rf_classification\257273\sampled_maps\';
for i_rf = 1:8
    
    for i_tree = 1:9
    
        sampled_pts_name = ['sampled_pts' zerostr(curr_tree, 4), '.mat'];
        sampled_pts_list{curr_tree} = [sampled_maps_dir sampled_pts_name]; %#ok

        copyfile(...
            ['C:\isbe\nailfold\models\vessel\detection\rf_classification\257273\01_sampled_data\sampled_pts' zerostr(i_tree,3) '.mat'],...
            sampled_pts_list{curr_tree});
        curr_tree = curr_tree + 1;
    end
end

make_sampled_maps(image_list, sampled_pts_list, sampled_maps_dir, 1);
%%
apex_radius = 5;
max_apex_guess = 500;
vessel_prob_smoothing_sigma = 2;
curvature_smoothing_sigma = 2;
strong_vessel_thresh = 0.25;
curv_max = 0.5;

g = gaussian_filters_1d(vessel_prob_smoothing_sigma);
g = g / sum(g);

rsa_dir = 'rsa_study/';
test_dir = 'set12';

apex_dir = [nailfoldroot 'data/' rsa_dir 'apexes/'];
model_dir = [nailfoldroot 'data/' rsa_dir 'models/apex_templates/'];
image_dir = [nailfoldroot 'data/' rsa_dir test_dir '/images/'];
fov_mask_dir = [nailfoldroot 'data/' rsa_dir test_dir '/fov_masks/'];
prob_dir = [nailfoldroot 'data/' rsa_dir test_dir '/predictions/detection/rf_classification/257273/'];
ori_dir = [nailfoldroot 'data/' rsa_dir test_dir '/predictions/orientation/rf_regression/259076/'];
width_dir = [nailfoldroot 'data/' rsa_dir test_dir '/predictions/width/rf_regression/257847/'];

pred_list = dir([prob_dir '*.mat']);
%%
patch_sz = 64;
patch_sz = patch_sz + 2; %Account for padding
patch_sz2 = (patch_sz - 1)/2;

x = repmat(-patch_sz2:patch_sz2, patch_sz, 1);
y = x';
xy = [x(:) y(:)];
dist_thresh = 24^2;
base_width = 20;

for i_im = 251:270
    
    vessel_im = u_load([image_dir pred_list(i_im).name(1:end-9) '.mat']);
    vessel_prob = u_load([prob_dir pred_list(i_im).name]);
    vessel_ori = u_load([ori_dir pred_list(i_im).name]);
    width_map = u_load([width_dir pred_list(i_im).name]);
    
    vessel_prob = conv2(g', g, vessel_prob, 'same');
    width_map = conv2(g', g, width_map, 'same');
    vessel_ori = conv2(g', g, vessel_ori, 'same');
        
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
    
    nailfold_vec = compute_weighted_histogram(width_map, vessel_prob, 0:2:150);
    
    if pred_list(i_im).name(1) == 'e'
        apex_struc = load([apex_dir 'enlarged/' pred_list(i_im).name(9:16) '.mat']);
        vessel_struc = u_load([apex_dir 'enlarged/' pred_list(i_im).name(9:16) '_vessel.mat']);
    elseif pred_list(i_im).name(1) == 'n'
        apex_struc = load([apex_dir 'normal/' pred_list(i_im).name(7:14) '.mat']);
        vessel_struc = u_load([apex_dir 'normal/' pred_list(i_im).name(7:14) '_vessel.mat']);
    elseif pred_list(i_im).name(1) == 'g'
        apex_struc = load([apex_dir 'giant/' pred_list(i_im).name(6:13) '.mat']);
        vessel_struc = u_load([apex_dir 'giant/' pred_list(i_im).name(6:13) '_vessel.mat']);
    end

    %Correct apex coordinates frame
    apex_xy(:,1) = apex_struc.apex_xy(:,1) +...
        apex_struc.apex_properties.sc - vessel_struc.vessel_properties.sc;
    apex_xy(:,2) = apex_struc.apex_xy(:,2) +...
        apex_struc.apex_properties.sr - vessel_struc.vessel_properties.sr;
    
    width_gt = sqrt(sum(diff(apex_xy).^2));
    
    dists = (vessel_centre_x - mean(apex_xy(:,1))).^2 +...
        (vessel_centre_y - mean(apex_xy(:,2))).^2;
    [~, min_idx] = min(dists);
    
    vxc = vessel_centre_x(min_idx);
    vyc = vessel_centre_y(min_idx);
    ori_c = angle(vessel_ori(vyc, vxc))/2;
    width_c = width_map(vyc, vxc);
    
    rot = [cos(ori_c) -sin(ori_c); sin(ori_c) cos(ori_c)];
    scale = width_c / base_width;

    %Transform points given scale and angle and translate to
    %candidate position
    xya = xy * rot * scale;
    xa = reshape(xya(:,1) + vxc, patch_sz, patch_sz);
    ya = reshape(xya(:,2) + vyc, patch_sz, patch_sz);
    
    vessel_prob_patch = interp2(vessel_prob, xa, ya, '*linear', 0);
    vessel_im_patch = interp2(vessel_im, xa, ya, '*linear', 0);
    
    [hog] = compute_HoG(vessel_prob_patch);
    
    figure;
    subplot(2,2,1); imgray(complex2rgb(vessel_ori));
    plot(xa(:,[1 end]), ya(:,[1 end]), 'b.');
    plot(xa([1 end],:), ya([1 end],:), 'b.');
    title('Vessel probability');
    
    subplot(2,2,3:4); bar(0:2:150, nailfold_vec(:));%imgray(vessel_prob);
    title('Histogram of widths weighted by vessel probability');
    
    plot_HoG(vessel_prob_patch, hog, [8 8], 1, subplot(2,2,2));
    title(['Predicted width: ' num2str(round(width_c)) ', annotated width: ' num2str(round(width_gt))]);
    
    %plot(vessel_centre_x, vessel_centre_y, 'b.');
    %plot(apex_xy(:,1), apex_xy(:,2), 'b');
    %plot(box_inner(:,1), box_inner(:,2), 'b');
    
    
    
end
%%
apex_radius = 5;
max_apex_guess = 500;
vessel_prob_smoothing_sigma = 2;
curvature_smoothing_sigma = 0.25;
strong_vessel_thresh = 0.25;
curv_max = 0.5;
smoothing_window = ones(1,25) / 25;

g = gaussian_filters_1d(vessel_prob_smoothing_sigma);
g = g / sum(g);

rsa_dir = 'rsa_study/';
test_dir = 'test';

apex_gt_dir = [nailfoldroot 'data/' rsa_dir test_dir '/apex_gt/'];
model_dir = [nailfoldroot 'data/' rsa_dir 'models/apex_templates/'];
image_dir = [nailfoldroot 'data/' rsa_dir test_dir '/images/'];
fov_mask_dir = [nailfoldroot 'data/' rsa_dir test_dir '/fov_masks/'];
prob_dir = [nailfoldroot 'data/' rsa_dir test_dir '/predictions/detection/rf_classification/257273/'];
ori_dir = [nailfoldroot 'data/' rsa_dir test_dir '/predictions/orientation/rf_regression/222835/'];
width_dir = [nailfoldroot 'data/' rsa_dir test_dir '/predictions/width/rf_regression/257847/'];

pred_list = dir([prob_dir '*.mat']);
%
num_images = length(pred_list);
%%
total_ims = 1;
%
max_ims = 20;
for i_im = 1:num_images
    
    im_num = pred_list(i_im).name(1:6);
    
    try
        vessel_prob = u_load([prob_dir pred_list(i_im).name]);
        vessel_ori = u_load([ori_dir pred_list(i_im).name]);
        vessel_width = u_load([width_dir pred_list(i_im).name]);
    catch
        continue;
    end
    vessel_im = u_load([image_dir im_num '.mat']);
    
    f_mask = u_load([fov_mask_dir im_num '_f_mask.mat']);
    load([apex_gt_dir im_num '_gt.mat'], 'vessel_markup');
    
    vessel_prob = imresize(vessel_prob, 0.5);
    vessel_ori = imresize(vessel_ori, 0.5);
    vessel_width = imresize(vessel_width, 0.5);
    f_mask = imresize(f_mask, 0.5);
    vessel_im = imresize(vessel_im,0.5);
    
    f_mask = imerode(f_mask, strel('disk', 50));
    
    vessel_prob(~f_mask) = 0;
    
    vessel_prob = conv2(g', g, vessel_prob, 'same');
    vessel_width = conv2(g', g, vessel_width, 'same');
    
    %Compute NMS centrelines
    vessel_nms = mb_non_maximal_supp(vessel_prob, angle(vessel_ori)/2);
    strong_vessels = vessel_nms > strong_vessel_thresh;
    if any(strong_vessels(:))
        [rstrong cstrong] = find(strong_vessels);
        vessel_centre_mask = bwselect(vessel_nms > 0, cstrong, rstrong, 8);
    else
        vessel_centre_mask = strong_vessels;
    end
    
    
    width_pdf = compute_weighted_histogram(vessel_width(vessel_centre_mask), vessel_prob(vessel_centre_mask), 0:10:150);
    
    vessel_prob_sm = conv2(smoothing_window', smoothing_window, vessel_prob, 'same');
    vessel_ori_sm = conv2(smoothing_window', smoothing_window, (vessel_ori ./ (abs(vessel_ori)+eps)).*vessel_prob_sm, 'same');
    vessel_ori_D = abs(vessel_ori_sm ./ vessel_prob_sm);
    vessel_ori_D(vessel_ori_D > 1) = 1;
    vessel_ori_D(isnan(vessel_ori_D)) = 1;
    ori_D_pdf = compute_weighted_histogram(vessel_ori_D(vessel_centre_mask), vessel_prob(vessel_centre_mask), 0:0.05:1);

    figure; 
    subplot(2,2,1:2); imgray(vessel_im); caxis([min(vessel_im(f_mask)) max(vessel_im(f_mask))]);
    title([im_num ': ' vessel_markup.image_grade ' graded by ' vessel_markup.observer]);
    
    subplot(2,2,3); bar(1:length(width_pdf(:)), width_pdf(:));
    subplot(2,2,4); bar(1:length(ori_D_pdf(:)), ori_D_pdf(:));
    
    total_ims = total_ims + 1;
    if total_ims > max_ims
        break;
    end
end
%%
