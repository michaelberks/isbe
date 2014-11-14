d_root = [asymmetryroot 'data/retinograms/dRIVE/test/'];
vessel_prob = load_uint8([d_root 'predictions/centre_detection/dt/rf_3/' zerostr(1,2) '_test_ext_class.mat']);

seed_map = zeros(size(vessel_prob));
for ii = 1:5e5
    %Choose a starting point - possibly provide options on this?
    [y_pot x_pot] = find(vessel_prob  > rand);
    r_idx = ceil(rand*length(x_pot));
    x1 = x_pot(r_idx);
    y1 = y_pot(r_idx);
    
    seed_map(y1, x1) = seed_map(y1, x1) + 1;
end
%%
for ii = 1:20
    path_map_co = u_load(['C:\isbe\asymmetry_project\data\retinograms\DRIVE\test\predictions\prob_paths\no_j\' zerostr(ii,2) '_path_map_c.mat']);
    path_map_ca = u_load(['C:\isbe\asymmetry_project\data\retinograms\DRIVE\test\predictions\prob_paths\no_j\' zerostr(ii,2) '_path_c_a.mat']);
    path_map_cb = u_load(['C:\isbe\asymmetry_project\data\retinograms\DRIVE\test\predictions\prob_paths\no_j\' zerostr(ii,2) '_path_c_a.mat']);
    
    figure; 
    
    path_map_c = path_map_cb + path_map_co;
    save(['C:\isbe\asymmetry_project\data\retinograms\DRIVE\test\predictions\prob_paths\no_j\10e5\' zerostr(ii,2) '_path_map_c.mat'], 'path_map_c');
    subplot(1,2,1); imgray(path_map_c);
    
    path_map_c = path_map_ca + path_map_cb + path_map_co;
    save(['C:\isbe\asymmetry_project\data\retinograms\DRIVE\test\predictions\prob_paths\no_j\12e5\' zerostr(ii,2) '_path_map_c.mat'], 'path_map_c');
    subplot(1,2,2); imgray(path_map_c);
    
%     path_map_c = path_map_ca + path_map_co;
%     save(['C:\isbe\asymmetry_project\data\retinograms\DRIVE\test\predictions\prob_paths\no_j\7.5e5\' zerostr(ii,2) '_path_map_c.mat'], 'path_map_c');
%     subplot(1,3,3); imgray(path_map_c);
end
%%
CC = bwconncomp(BW);
numPixels = cellfun(@numel,CC.PixelIdxList);
[biggest,idx] = max(numPixels);
BW(CC.PixelIdxList{idx}) = 0;
%%
jj = 1;    
data = 'test';
win_size = 11;
    
%Load in ground truth and skeletonise
gt = imread(['C:\isbe\asymmetry_project\data\retinograms\DRIVE\' data '\1st_manual\' zerostr(jj,2) '_manual1.gif']);
gt = gt > 0;
gts = bwmorph(gt, 'skel', 'inf');


%     figure; 
%     subplot(2,2,1); imagesc(gt); axis image;
%     subplot(2,2,2); imagesc(gts); axis image;

%Extract x,y coords of vessels + vessel centres
[c_y c_x] = find(gts);
[a_y a_x] = find(gt);

%Create storage for the ground truth orientations
gts_ori = zeros(size(gts));

%Loop through each skeleton point
for ii = 1:size(c_x,1)

    %Sample local window from skeleton map
    local_win = sample_window(gts, win_size, c_y(ii), c_x(ii), 0);

    %Get all points connected to the centre
    [yi xi] = find(bwselect(local_win, (win_size+1)/2, (win_size+1)/2, 8));
    uni_x = unique(xi);
    uni_y = unique(yi);

    if length(uni_x) > length(uni_y)
        uni_y = sparse(xi, 1, yi, win_size, 1) ./ sparse(xi, 1, 1, win_size, 1);
        uni_y = full(uni_y(uni_x));
    else
        uni_x = sparse(yi, 1, xi, win_size, 1) ./ sparse(yi, 1, 1, win_size, 1);
        uni_x = full(uni_x(uni_y));
    end

    uu = mean(diff(uni_x));
    vv = -mean(diff(uni_y));
    dd = sqrt(uu^2 + vv^2);
    gts_ori(c_y(ii), c_x(ii)) = complex(uu / dd, vv / dd);
end
%%
a_u = griddata(c_x, c_y, real(gts_ori(gts)),a_x, a_y, 'nearest');
a_v = griddata(c_x, c_y, imag(gts_ori(gts)),a_x, a_y, 'nearest');

dd = a_u.^2 + a_v.^2;

gt_ori = zeros(size(gt));
gt_ori(gt) = (complex(a_u, a_v).^2) ./ dd;
%%
r = 300:334;
c = 357:391;
gts_roi = gts(r, c);
ori_roi = gts_ori(r, c);

[y x] = find(gts_roi);
theta = ori_roi(gts_roi);
figure; imgray(gts_roi);
quiver(x,y,2*real(theta),-2*imag(theta),'ShowArrowHead', 'off','Autoscale', 'off');
quiver(x,y,-2*real(theta),2*imag(theta),'ShowArrowHead', 'off','Autoscale', 'off');
%%
%now loop through gts again...
num_pts = size(c_x,1);
gts_ori_mixed = cell(num_pts,1);
win_size_m = 7;
gt_ori_idx = zeros(size(gts));
for ii = 1:num_pts

    %Sample local window from skeleton map
    local_win = sample_window(gts, win_size_m, c_y(ii), c_x(ii), 0);
    local_ori = sample_window(gt_ori, win_size_m, c_y(ii), c_x(ii), 0);
    local_connections = bwselect(local_win, (win_size_m+1)/2, (win_size_m+1)/2, 8);

    gts_ori_mixed{ii} = local_ori(local_connections);
end
a_idx = griddata(c_x, c_y, 1:num_pts,a_x, a_y, 'nearest');
gt_ori_idx(gt) = a_idx;

%%
colors = lines(12);
figure; imgray(gt);
for ii = 1:size(a_x,1)
    x = a_x(ii);
    y = a_y(ii);
    
    if x < 340 || x > 400 || y < 320 || y > 380
        continue;
    end
    idx = gt_ori_idx(y,x);
    num_oris = length(gts_ori_mixed{idx});
    mixed_oris = gts_ori_mixed{idx};
    mixed_oris = abs(mixed_oris) .* exp(1i*angle(mixed_oris)/2); 
    color = mod(idx,12) + 1;
    quiver(ones(num_oris,1)*x,ones(num_oris,1)*y,real(mixed_oris),-imag(mixed_oris), 'color', colors(color,:));
    quiver(ones(num_oris,1)*x,ones(num_oris,1)*y,-real(mixed_oris),imag(mixed_oris), 'color', colors(color,:));
end
%%
[training_data training_labels] = generate_vessel_data( ... % non-strict mode
    'num_samples', 1e4,...
    'image_dir', 'C:\isbe\asymmetry_project\data\retinograms\DRIVE\training\images_extended\',...
    'foveal_mask_dir','C:\isbe\asymmetry_project\data\retinograms\DRIVE\training\foveal_masks\',...
    'vessel_mask_dir','C:\isbe\asymmetry_project\data\retinograms\DRIVE\training\vessel_masks\',... % the mandatory arguments
    'prediction_type', 'mixed_orientation',...
    'ori_dir', 'C:\isbe\asymmetry_project\data\retinograms\DRIVE\training\mixed_orientations\',...
    'selected_images', 1,...
    'rgb_channel', 'rgb',...
    'win_size', 3,...
    'num_levels', 4,...
    'feature_type', 'conj',...
    'feature_shape', 'rect',...
    'decomp_type', 'dt',...
    'bg_ratio', 1);
%%
[training_data training_labels] = generate_vessel_data( ... % non-strict mode
    'num_samples', 5e3,...
    'image_dir', 'C:\isbe\asymmetry_project\data\retinograms\DRIVE\training\images\',...
    'foveal_mask_dir','C:\isbe\asymmetry_project\data\retinograms\DRIVE\training\foveal_masks\',...
    'vessel_mask_dir','C:\isbe\asymmetry_project\data\retinograms\DRIVE\training\vessel_masks\',... % the mandatory arguments
    'prediction_type', 'foveal_detection',...
    'selected_images', 1,...
    'rgb_channel', 'rgb',...
    'win_size', 1,...
    'decomp_type', 'pixel',...
    'bg_ratio', 1);

%%
d_root = [asymmetryroot 'data/retinograms/dRIVE/test/'];
ret = rgb2gray(imread([d_root 'images\01_test.tif']));
v_mask = u_load([d_root 'vessel_masks\' zerostr(1,2) '_test_v_mask.mat']);
vs_mask = bwmorph(v_mask, 'skel', 'inf');
figure; imgray(vs_mask);

%%
forest = u_load('C:\isbe\asymmetry_project\data\models\vessel\centre_orientation\73939\random_forest');
forest.tree_root = 'C:\isbe\asymmetry_project\data\models\vessel\centre_orientation\';
sampling_args = u_load('C:\isbe\asymmetry_project\data\models\vessel\centre_orientation\73939\sampling_args.mat');
sampling_args_c.num_levels = sampling_args.num_levels;
sampling_args_c.feature_shape = sampling_args.feature_shape;
sampling_args_c.feature_type = sampling_args.feature_type;
sampling_args_c.do_max = sampling_args.do_max;
sampling_args_c.rotate = sampling_args.rotate;
sampling_args_c.win_size = sampling_args.win_size;
if isfield(sampling_args, 'use_nag')
    sampling_args_c.use_nag = sampling_args.use_nag;
else
    sampling_args_c.use_nag = args.use_nag;
end
%%
forest = u_load('C:\isbe\asymmetry_project\data\models\vessel\mixed_orientation\77098\random_forest');
forest.tree_root = 'C:\isbe\asymmetry_project\data\models\vessel\mixed_orientation\';
sampling_args = u_load('C:\isbe\asymmetry_project\data\models\vessel\mixed_orientation\77098\sampling_args.mat');
sampling_args_c.num_levels = sampling_args.num_levels;
sampling_args_c.feature_shape = sampling_args.feature_shape;
sampling_args_c.feature_type = sampling_args.feature_type;
sampling_args_c.do_max = sampling_args.do_max;
sampling_args_c.rotate = sampling_args.rotate;
sampling_args_c.win_size = sampling_args.win_size;
if isfield(sampling_args, 'use_nag')
    sampling_args_c.use_nag = sampling_args.use_nag;
else
    sampling_args_c.use_nag = args.use_nag;
end
%%
j_mask = false(size(v_mask));
j_mask(272,313) = 1;
classify_image(...
    'image_in',ret,... % the mandatory arguments
    'sampling_args',sampling_args_c,...
    'forest', forest, ...
    'decomp_type', 'dt',...
    'forest_type', 'orientation',...
    'use_probs', 0,...
    'mask', j_mask,...
    'num_trees', [], ...
    'max_size', 128);
%%
uv = abs(y_trees).*exp(1i*angle(y_trees)/2);
plot([313 - 5*real(uv); 313 + 5*real(uv)], [ 272 + 5*imag(uv); 272 - 5*imag(uv)]);
plot([313 - 5*real(mean(uv)); 313 + 5*real(mean(uv))], [ 272 + 5*imag(mean(uv)); 272 - 5*imag(mean(uv))], 'g', 'linewidth', 4);
%%
d_root = [asymmetryroot 'data/retinograms/dRIVE/test/'];

centre_vessels = 0;
nms = 0;

for roc_method = [1]
    
    tp_counts = zeros(20, 102);
    fp_counts = zeros(20, 102);
    t_counts = zeros(20, 1);
    f_counts = zeros(20, 1);

    f1 = figure; axis equal; axis([0 1 0 1]); hold all;
    xlabel('FPF');
    ylabel('TPF');
    for data_type = 1
        for ii = 1:20

            switch data_type

                case 1
                    label = 'Raw classifications, centre-line, non-dilated';
                    vessel_prob = load_uint8(['Z:\asym\data\retinograms\DRIVE\test\predictions\detection\dt\rf_3/' zerostr(ii,2) '_test_class.mat']);

                case 2
                    label = 'Staal';
                    vessel_prob = double(rgb2gray(imread([d_root 'predictions\staal\' zerostr(ii-1, 2) '.bmp'])))/255; 
%                     label = 'Tracking paths, centre-line';
%                     load([d_root 'predictions/prob_paths/good/' zerostr(ii,2) '_path_c.mat'], 'path_map_c');
%                     vessel_prob = log(path_map_c) ./ max(log(path_map_c(:)));%prctile(path_map_c(:), 99.5);
%                     vessel_prob(vessel_prob == -inf) = 0;
                case 3
                    label = 'Tracking paths, centre-line';
                    load([d_root 'predictions/prob_paths/thin/' zerostr(ii,2) '_path_c.mat'], 'path_map_c');
                    vessel_prob = log(path_map_c) ./ max(log(path_map_c(:)));%prctile(path_map_c(:), 99.5);
                    vessel_prob(vessel_prob == -inf) = 0;
                case 4
                    label = 'Raw classifications, centre-line';
                    vessel_prob = load_uint8([d_root 'predictions/centre_detection/dt/rf_3/' zerostr(ii,2) '_test_ext_class.mat']);
                    
            end

            if nms
                [~, vessel_theta] = gaussian_2nd_derivative_line(vessel_prob, 1.5);
                vessel_prob = mb_non_maximal_supp(vessel_prob, vessel_theta);
            end

            v_mask = u_load([d_root 'vessel_masks\' zerostr(ii,2) '_test_v_mask.mat']);
            f_mask = u_load([d_root 'foveal_masks\' zerostr(ii,2) '_test_f_mask.mat']);
            vs_mask = v_mask;%bwmorph(v_mask, 'thin', 'inf');
            
            %Compute ROC counts for image
            switch roc_method
                case 1
                    [~, ~, tp_count fp_count] = calculate_roc_image(vessel_prob, vs_mask,(-1:100)/100, f_mask, 'dilate', 0);
                    t_counts(ii) = sum(vs_mask(f_mask));
                    f_counts(ii) = sum(~vs_mask(f_mask));
                    
                case 2
                    [~, ~, tp_count fp_count] = calculate_roc_image(vessel_prob, vs_mask,(-1:100)/100, f_mask, 'dilate', 1);
                    t_counts(ii) = sum(vs_mask(f_mask));
                    f_counts(ii) = sum(~vs_mask(f_mask));
                    
                case 3
                    [~, ~, tp_count fp_count] = calculate_roc_image2(vessel_prob, vs_mask & f_mask, ~v_mask & f_mask, (-1:100)/100, 'dilate', 0);
                    t_counts(ii) = sum(vs_mask(f_mask));
                    f_counts(ii) = sum(~v_mask(f_mask));
                case 4
                    [~, ~, tp_count fp_count] = calculate_roc_image2(vessel_prob, vs_mask & f_mask, ~v_mask & f_mask, (-1:100)/100, 'dilate', 1);
                    t_counts(ii) = sum(vs_mask(f_mask));
                    f_counts(ii) = sum(~v_mask(f_mask));
                case 5
                    [~, ~, tp_count fp_count] = calculate_roc_connected(vessel_prob, vs_mask,(-1:100)/100, f_mask);
                    t_counts(ii) = sum(vs_mask(f_mask));
                    f_counts(ii) = sum(~vs_mask(f_mask));
                
            end

            %Increment total counts
            tp_counts(ii,:) = tp_count;
            fp_counts(ii,:) = fp_count;
                

        end

        %Compute ROC points for complete set of data
        roc_pts = [sum(fp_counts)' / sum(f_counts) sum(tp_counts)' / sum(t_counts)];

        %Compute AUC for ROC curve
        auc = sum( (roc_pts(2:end,1)-roc_pts(1:end-1,1)) .* roc_pts(1:end-1,2)) + ...
                0.5*sum( (roc_pts(2:end,1)-roc_pts(1:end-1,1)) .* (roc_pts(2:end,2)-roc_pts(1:end-1,2)) );


        figure(f1);
        plot(roc_pts(:,1), roc_pts(:,2), '-'); axis([0 1 0 1]);
        plot(roc_pts(:,1), roc_pts(:,2), '.');
        legend_label{2*data_type-1,1} = ' ';
        legend_label{2*data_type,1} = [ label ' AUC: ' num2str(auc)]; %#ok
    end
    legend(legend_label, 'location', 'southeast');
end
%%
%Accuracy table
accuracy = zeros(20,100);
for ii = 1:20
    for jj = 1:100
        vessel_prob = load_uint8([d_root 'test\images/results/77469/' zerostr(ii,2) '_test_class.mat']);
        v_mask = u_load([d_root 'test\vessel_masks\' zerostr(ii,2) '_test_v_mask.mat']);
        f_mask = u_load([d_root 'test\foveal_masks\' zerostr(ii,2) '_test_f_mask.mat']);
        hard_class = vessel_prob > (jj/100);
        accuracy(ii,jj) = sum(hard_class(f_mask) == v_mask(f_mask)) / sum(f_mask(:));
    end
end
mean_accuracies = mean(accuracy);
[max_accuracy operating_pt] = max(mean_accuracies)
%%
accuracy = zeros(20,100);
for ii = 1:20
    for jj = 1:100
        load([d_root 'predictions/prob_paths/no_edge/' zerostr(ii,2) '_path.mat'], 'path_map');
        vessel_prob = log(path_map) ./ max(log(path_map_c(:)));%prctile(path_map_c(:), 99.5);
        vessel_prob(vessel_prob == -inf) = 0;
        
        v_mask = u_load([d_root 'vessel_masks\' zerostr(ii,2) '_test_v_mask.mat']);
        f_mask = u_load([d_root 'foveal_masks\' zerostr(ii,2) '_test_f_mask.mat']);
        hard_class = vessel_prob > (jj/100);
        accuracy(ii,jj) = sum(hard_class(f_mask) == v_mask(f_mask)) / sum(f_mask(:));
    end
end
mean_accuracies = mean(accuracy);
[max_accuracy operating_pt] = max(mean_accuracies)
%%
d_root = [asymmetryroot 'data/retinograms/dRIVE/test/'];
for ii = 1:20
     ret = imread([d_root 'images\' zerostr(ii,2) '_test.tif']);
     save([d_root 'images\' zerostr(ii,2) '_test.mat'], 'ret');
end
d_root = [asymmetryroot 'data/retinograms/dRIVE/training/'];
for ii = 21:40
     ret = imread([d_root 'images\' zerostr(ii,2) '_training.tif']);
     save([d_root 'images\' zerostr(ii,2) '_training.mat'], 'ret');
end
%%
clear
d_root = [asymmetryroot 'data/retinograms/dRIVE/'];

forest = u_load([asymmetryroot 'data\models\vessel\foveal_detection\77325\random_forest']);
forest.tree_root = [asymmetryroot 'data\models\vessel\foveal_detection\'];
sampling_args = u_load([asymmetryroot '\models\vessel\foveal_detection\77325\sampling_args.mat']);
sampling_args_c.num_levels = sampling_args.num_levels;
sampling_args_c.feature_shape = sampling_args.feature_shape;
sampling_args_c.feature_type = sampling_args.feature_type;
sampling_args_c.do_max = sampling_args.do_max;
sampling_args_c.rotate = sampling_args.rotate;
sampling_args_c.win_size = sampling_args.win_size;
if isfield(sampling_args, 'use_nag')
    sampling_args_c.use_nag = sampling_args.use_nag;
else
    sampling_args_c.use_nag = args.use_nag;
end

ret_te1 = u_load([d_root 'test/temp/01_test.mat']);
f_mask = u_load([d_root 'test/foveal_masks/01_test_f_mask.mat']);

classify_image(...
    'image_in',rgb2gray(ret_te1),... % the mandatory arguments
    'sampling_args',sampling_args_c,...
    'forest', forest, ...
    'decomp_type', 'dt',...
    'forest_type', 'isbe',...
    'use_probs', 0,...
    'mask', f_mask,...
    'num_trees', [], ...
    'max_size', 128);
%%
clear
d_root = [asymmetryroot 'data/retinograms/dRIVE/'];

forest = u_load([asymmetryroot 'data\models\vessel\detection\51101\random_forest']);
forest.tree_root = [asymmetryroot 'data\models\vessel\detection\'];
sampling_args = u_load([asymmetryroot 'data\models\vessel\detection\51101\sampling_args.mat']);
sampling_args_c.num_levels = sampling_args.num_levels;
sampling_args_c.feature_shape = sampling_args.feature_shape;
sampling_args_c.feature_type = sampling_args.feature_type;
sampling_args_c.do_max = sampling_args.do_max;
sampling_args_c.rotate = sampling_args.rotate;
sampling_args_c.win_size = sampling_args.win_size;
if isfield(sampling_args, 'use_nag')
    sampling_args_c.use_nag = sampling_args.use_nag;
else
    sampling_args_c.use_nag = args.use_nag;
end

ret_te1 = u_load([d_root 'test/temp/01_test.mat']);
f_mask = u_load([d_root 'test/foveal_masks/01_test_f_mask.mat']);

vp_goat = classify_image(...
    'image_in',rgb2gray(ret_te1),... % the mandatory arguments
    'sampling_args',sampling_args_c,...
    'forest', forest, ...
    'decomp_type', 'dt',...
    'forest_type', 'isbe',...
    'use_probs', 0,...
    'mask', f_mask,...
    'num_trees', [], ...
    'max_size', 128);
%%
%% 1. Produce response, orientation and scale maps for the analytic methods
retroot = [asymmetryroot('shared'),'data\retinograms\DRIVE\test\predictions\'];
do_mono =0;
do_g1d = 0;
do_g2d = 1;

if do_mono
    mkdir([retroot 'orientations\mono\analytic\orientations']);
    mkdir([retroot 'orientations\mono\analytic\responses']);
    mkdir([retroot 'orientations\mono\analytic\scales']);
end
if do_g1d
    mkdir([retroot 'orientations\g1d\analytic\orientations']);
    mkdir([retroot 'orientations\g1d\analytic\responses']);
    mkdir([retroot 'orientations\g1d\analytic\scales']);
end
if do_g2d
    mkdir([retroot 'orientations\g2d_o\analytic\orientations']);
    mkdir([retroot 'orientations\g2d_o\analytic\responses']);
    mkdir([retroot 'orientations\g2d_o\analytic\scales']);
end

for ii = 1:20
    %load retinogram and merge RGB channels
    ret = u_load([asymmetryroot('shared'),'data\retinograms\DRIVE\test\images\' zerostr(ii,2) '_test.mat']);
    ret = rgb2gray(ret);
    
    if do_mono
        %Compute repsonses for mono and save
        [response_map d ori_map scale_map] = monogenic_multiscale(ret, 4, 4, 2, 0.65);
        save_uint8([retroot 'mono\analytic\orientations\' zerostr(ii,3) '_test_ori.mat'], ori_map);
        save_uint8([retroot 'mono\analytic\responses\' zerostr(ii,3) '_test_response.mat'], response_map);
        save_uint8([retroot 'mono\analytic\scales\' zerostr(ii,3) '_test_scale.mat'], scale_map);
    end
    if do_g1d
        %Compute repsonses for mono and save
        [response_map ori_map scale_map] = gaussian_1st_derivative_gradient2(ret, [1 2 4 8]);
        save_uint8([retroot 'g1d\analytic\orientations\' zerostr(ii,3) '_test_ori.mat'], ori_map);
        save_uint8([retroot 'g1d\analytic\responses\' zerostr(ii,3) '_test_response.mat'], response_map);
        save_uint8([retroot 'g1d\analytic\scales\' zerostr(ii,3) '_test_scale.mat'], scale_map);
    end
    if do_g2d
        %Compute repsonses for mono and save
        [response_map ori_map scale_map] = gaussian_clover_line(ret, [1 2 4 8]);
        save_uint8([retroot 'orientations\g2d_o\analytic\orientations\' zerostr(ii,3) '_test_ori.mat'], ori_map);
        save_uint8([retroot 'orientations\g2d_o\analytic\responses\' zerostr(ii,3) '_test_response.mat'], response_map);
        save_uint8([retroot 'orientations\g2d_o\analytic\scales\' zerostr(ii,3) '_test_scale.mat'], scale_map);
    end
    

end
%%

mkdir([retroot 'orientations\post_classc_g2d\analytic\orientations']);
mkdir([retroot 'orientations\post_classc_g2d\analytic\responses']);
mkdir([retroot 'orientations\post_classc_g2d\analytic\scales']);

for ii = 1:20
    %load retinogram and merge RGB channels
    ret = load_uint8([asymmetryroot 'data\retinograms\DRIVE\test\images/results/77541/' zerostr(ii,2) '_test_class.mat']);

    %Compute repsonses for mono and save
    [response_map ori_map scale_map] = gaussian_clover_line(ret, [1 2 4 8]);
    save_uint8([retroot 'orientations\post_classc_g2d\analytic\orientations\' zerostr(ii,3) '_test_ori.mat'], ori_map);
    save_uint8([retroot 'orientations\post_classc_g2d\analytic\responses\' zerostr(ii,3) '_test_response.mat'], response_map);
    save_uint8([retroot 'orientations\post_classc_g2d\analytic\scales\' zerostr(ii,3) '_test_scale.mat'], scale_map);

end

%% 2. Now compute error statistics for all the methods
retroot = [asymmetryroot,'data\retinograms\DRIVE\test\'];

%first make sure we've got a complete set of ground truth labels
if ~exist([retroot,'orientations\all_gt_orientations.mat'],'file')
	gt_orientations = [];
	vessel_centres = [];
	pts_per_image = zeros(20,1); 
	for jj = 1:20

		%load foveal mask and vessel mask and indices to vessel centres
		foveal_mask = u_load([retroot 'foveal_masks\' zerostr(jj,2) '_test_f_mask']);
		vessel_mask = u_load([retroot 'vessel_masks\' zerostr(jj,2) '_test_v_mask.mat']);
		vessel_mask(~foveal_mask) = false;

		vessel_centres = logical([vessel_centres; ...
			u_load([retroot 'vessel_masks\centre_idx\' zerostr(jj,2) '_test_centre_idx.mat'])]);

		%load in the ground truth orientation
		gt_ori = u_load([retroot,'orientations\' zerostr(jj,2) '_ori1.mat']);
		gt_orientations = [gt_orientations; gt_ori(vessel_mask)]; %#ok
		pts_per_image(jj) = sum(vessel_mask(:));
	end
	save([retroot,'orientations\all_gt_orientations.mat'], 'gt_orientations', 'pts_per_image', 'vessel_centres');
end
%%
