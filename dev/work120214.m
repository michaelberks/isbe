jetstream(vessel_prob, gd1_or, [60 6], [cos(gd1_or(6,60)) -sin(gd1_or(6,60))],...
    'sigma_theta', 0.1, 'N', 200, 'plot', 1);
%%
ii = 1;
d_root = 'C:\isbe\asymmetry_project\data\retinograms\dRIVE\test\';
%vessel_prob = load_uint8([d_root 'images_extended\results\73563\' zerostr(ii,2) '_test_ext_class.mat']);
vessel_prob = load_uint8([d_root 'images_extended\results\51101\' zerostr(ii,2) '_test_ext_class.mat']);
vessel_ori = load_uint8(['Z:\asym\data\retinograms\DRIVE\test\predictions2\dt\rf_3\' zerostr(ii,2) '_ori.mat']);
j_map = load_uint8(['Z:\asym\data\retinograms\DRIVE\test\images_extended\results\73722\198\' zerostr(ii,2) '_test_ext_class.mat']);
f_mask = u_load([d_root 'foveal_masks\' zerostr(ii,2) '_test_f_mask.mat']);
vessel_prob(~f_mask) = 0;

vessel_theta = angle(vessel_ori);
vessel_mag = abs(vessel_ori);
%%
%figure; imgray(vessel_prob);
path_map = zeros(size(vessel_prob));
for jj = 1:1e5
    xi = 321;% 81;
    yi = 102;%252;
    t1 = vessel_theta(round(yi), round(xi));%rand*pi - pi/2 ;
    xs0 = cos(t1);
    ys0 = -sin(t1);

    x = xi;
    y = yi;
    go_on = true;
    step = 1;
    while go_on

        is_j = j_map(round(yi), round(xi));
        if 0%is_j > .5
            new_theta = 2*pi*rand;
        else
            theta = vessel_theta(round(yi), round(xi));
            rho = vessel_mag(round(yi), round(xi));
            new_theta = theta + wrapped_normal_sample(0, rho, 1)/2;
        end
            
        xs1 = cos(new_theta);
        ys1 = -sin(new_theta);

        if xs1*xs0 + ys1*ys0 < -0.5
            xs1 = -xs1;
            ys1 = -ys1;
        end

        xi = xi + 2*xs1;
        yi = yi + 2*ys1;
        xs0 = xs1;
        ys0 = ys1;

        %x = [x; xi];
        %y = [y; yi];
        path_map(round(yi), round(xi)) = path_map(round(yi), round(xi)) + 2^((step-1)/2);

        p = vessel_prob(round(yi), round(xi));
        go_on = (p > 0.5) || (p > rand);
        step = step + 1;
    end
    %plot(x, y);
end
figure; imagesc(log(path_map)); axis image; colorbar;
% subplot(1,2,1); imagesc(path_map); axis image; colorbar;
% subplot(1,2,2); 
%%
total_pts = 0;

for jj = 1:20
    gt = u_load(['C:\isbe\asymmetry_project\data\retinograms\DRIVE\training\vessel_masks\' zerostr(20+jj,2) '_training_v_mask.mat']);
    %gts = bwmorph(gt, 'skel', 'inf');
    gtt = bwmorph(gt, 'thin', 'inf');
    v_thick = bwdist(~gt);
    [rows cols] = size(gt);
    xx = repmat(1:cols, rows, 1);
    yy = repmat((1:rows)', 1, cols);
    
    junctions = false(size(gt));
    
    [c_y c_x] = find(gtt);
    [a_y a_x] = find(gt);

    num_pts = length(c_x);

    ws = 5;
    ws2 = floor(ws/2);
    offs = -ws2:ws2;
    count_mask = ones(ws); count_mask(2:ws-1,2:ws-1) = 0;

    % offs = -1:1;
    % count_mask = ones(3); count_mask(2,2) = 0;
    is_branch = false(num_pts,1);
    %
    gtt2 = gtt;
    for ii = 1:num_pts

        patch = gtt(c_y(ii) + offs, c_x(ii) + offs);
        patch = bwselect(patch, ws2+1, ws2+1, 8) & count_mask;
        label = bwlabel(patch, 8);
        is_branch(ii) = max(label(:)) >= 3;
        if is_branch(ii)
            rad = v_thick(c_y(ii), c_x(ii));
            circ = (xx - c_x(ii)).^2 + (yy - c_y(ii)).^2 < rad^2;
            circ = circ & gt;
            junctions = junctions | circ;
        end
            
        gtt2(c_y(ii) + offs, c_x(ii) + offs) = max(label(:)) >= 3;

    end
    total_pts = total_pts + sum(gtt2(:));
%     figure; 
%     subplot(1,2,1); imgray(gt);
%     subplot(1,2,2); imgray(junctions);
    %plot(c_x(is_branch), c_y(is_branch), 'rx');
    
%     a_j = griddata(c_x, c_y, double(is_branch), a_x, a_y, 'nearest') > 0;
%     
%     junctions = double(gt);
%     junctions(gt) = junctions(gt) + a_j;
% 
%     figure; 
%     subplot(1,2,1); imagesc(gtt); axis image; colormap gray; hold on;
%     plot(c_x(is_branch), c_y(is_branch), 'rx');
% 
%     subplot(1,2,2); imagesc(junctions); axis image; colormap gray;
end
%%
total_pts = 0;
for jj = 1:20
    gt = u_load(['C:\isbe\asymmetry_project\data\retinograms\DRIVE\training\vessel_masks\' zerostr(20+jj,2) '_training_v_mask.mat']);
    %gts = bwmorph(gt, 'skel', 'inf');
    gtt = bwmorph(gt, 'thin', 'inf');
    
    total_pts = total_pts + sum(gtt(:));
end
%%
[training_data training_labels] = generate_vessel_data(...
    'num_samples', 1e5, ...
    'image_dir', 'Z:\asym\data\retinograms\DRIVE\training\images_extended\',...
    'foveal_mask_dir', 'Z:\asym\data\retinograms\DRIVE\training\foveal_masks\',...
    'vessel_mask_dir', 'Z:\asym\data\retinograms\DRIVE\training\vessel_masks\',... % the mandatory arguments
    'prediction_type', 'centre_orientation',...
    'ori_dir', 'Z:\asym\data\retinograms\DRIVE\training\orientations\',...
    'selected_images', [], ...
    'rgb_channel', 'g',...
    'win_size', 1,...
    'num_levels', 1,...
    'feature_type', 'conj',...
    'feature_shape', 'rect',...
    'rotate', 0,...
    'do_max', 0,...
    'decomp_type', 'pixel',...
    'bg_ratio', 1,...
    'use_nag', 0);
%%

for ii = 1:20
    v_mask = u_load(['Z:\asym\data\retinograms\DRIVE\test\vessel_masks\' zerostr(ii,2) '_test_v_mask.mat']);
    j_map1 = load_uint8(['Z:\asym\data\retinograms\DRIVE\test\images_extended\results\73722\198\' zerostr(ii,2) '_test_ext_class.mat']);
    j_map2 = load_uint8(['Z:\asym\data\retinograms\DRIVE\test\images_extended\results\73965\' zerostr(ii,2) '_test_ext_class.mat']);
    
    j_map_v1 = j_map1;
    j_map_v1(~v_mask) = 0;
    [jy1 jx1] = find(j_map_v1 > .5);
    
    j_map_v2 = j_map2;
    j_map_v2(~v_mask) = 0;
    [jy2 jx2] = find(j_map_v2 > .7);
    
    figure; 
    %subplot(1,2,1); imgray(j_map_v1); colormap hot;
    %subplot(1,2,2); imgray(j_map_v2); colormap hot;
    subplot(1,2,1); imgray(v_mask); plot(jx1, jy1, 'r.');
    subplot(1,2,2); imgray(v_mask); plot(jx2, jy2, 'r.');
end
%%
ii = 1;
d_root = 'C:\isbe\asymmetry_project\data\retinograms\dRIVE\test\';

%**** All vessel ******
vessel_prob = load_uint8([d_root 'predictions\detection\dt\rf_3\' zerostr(ii,2) '_test_ext_class.mat']);
vessel_ori = load_uint8([d_root 'predictions\orientation\dt\rf_3\' zerostr(ii,2) '_ori.mat']);

%**** Vessel centre ******
vessel_prob_c = load_uint8([d_root 'predictions\centre_detection\dt\rf_3\' zerostr(ii,2) '_test_ext_class.mat']);
vessel_ori_c = load_uint8([d_root 'predictions\centre_orientation\dt\rf_3\' zerostr(ii,2) '_test_ext_class.mat']);

j_map = load_uint8([d_root 'predictions\junction_detection\dt\rf_3\' zerostr(ii,2) '_test_ext_class.mat']);
j_map_c = load_uint8([d_root 'predictions\junction_centre_detection\dt\rf_1\' zerostr(ii,2) '_test_ext_class.mat']);

f_mask = u_load([d_root 'foveal_masks\' zerostr(ii,2) '_test_f_mask.mat']);
vessel_prob(~f_mask) = 0;
vessel_prob_c(~f_mask) = 0;

xi = 88;
yi = 245;
%[particles_x particles_y] = jetstream_rf(vessel_prob, vessel_ori, j_map>0.5, [xi yi], [sqrt(2) -sqrt(2)], [], 'plot', 1, 'N', 1000);
% [particles_x particles_y] = jetstream_rf(vessel_prob, vessel_ori, [], [xi yi], [sqrt(2) -sqrt(2)], [],...
%     'plot', 1, 'N', 1000, 'sigma_theta', 0.1);
[particles_x particles_y] = jetstream_rf2(vessel_prob, vessel_ori, j_map>1, [xi yi], [sqrt(2) -sqrt(2)], [],...
    'plot', 1, 'N', 1000, 'sigma_theta', 0.1);

%%
[path_map] = prob_track_mult(vessel_prob, vessel_ori, ...
    'num_paths', 1e5,...
    'junction_map', [],...
    'step_length', 2,...
    'double_angle', 0,...
    'plot', 1);

[path_map_c] = prob_track_mult(vessel_prob_c, vessel_ori_c, ...
    'num_paths', 1e5,...
    'junction_map', j_map_c,...
    'step_length', 2,...
    'double_angle', 1,...
    'plot', 1);
%%
for ii = 1%:20
    d_root = 'C:\isbe\asymmetry_project\data\retinograms\dRIVE\test\';

    %**** Vessel centre ******
    vessel_prob_c = load_uint8([d_root 'predictions\centre_detection\dt\rf_3\' zerostr(ii,2) '_test_ext_class.mat']);
    vessel_ori_c = load_uint8([d_root 'predictions\centre_orientation\dt\rf_3\' zerostr(ii,2) '_test_ext_class.mat']);

    f_mask = u_load([d_root 'foveal_masks\' zerostr(ii,2) '_test_f_mask.mat']);
    vessel_prob_c(~f_mask) = 0;
    
    group_map = staal_ridge_group(vessel_prob_c, vessel_ori_c);
    
    save([d_root 'predictions\ridges\' zerostr(ii,2) '_ridge.mat'], 'group_map');    
    
end
%%
for ii = 1%:20
    load([d_root 'predictions\ridges\' zerostr(ii,2) '_ridge.mat'], 'group_map');
    f_mask = u_load([d_root 'foveal_masks\' zerostr(ii,2) '_test_f_mask.mat']);
    
    [py px] = find(f_mask & ~group_map);
    [ry rx] = find(group_map > 0);
    
    set_nums = griddata(rx, ry, group_map(group_map > 0), px, py, 'nearest');
    
    set_map = group_map;
    set_map(f_mask & ~group_map) = set_nums;
    figure; imgray(set_map); colormap lines;
    
end
%%
d_root = [asymmetryroot 'data/retinograms/dRIVE/test/'];
mkdir([d_root 'predictions/prob_paths/']);
for ii = 1%:20
    

    %**** Vessel centre ******
    vessel_prob_c = load_uint8([d_root 'predictions/centre_detection/dt/rf_3/' zerostr(ii,2) '_test_ext_class.mat']);
    vessel_ori_c = load_uint8([d_root 'predictions/centre_orientation/dt/rf_3/' zerostr(ii,2) '_test_ext_class.mat']);

    f_mask = u_load([d_root 'foveal_masks/' zerostr(ii,2) '_test_f_mask.mat']);
    vessel_prob_c(~f_mask) = 0;
    
    [path_map_c] = prob_track_mult(vessel_prob_c, vessel_ori_c, ...
        'num_paths', 1e5,...
        'junction_map', [],...
        'ignore_mask', imerode(f_mask, strel('disk', 32)),...
        'step_length', 1,...
        'double_angle', 1,...
        'plot', 1);
    figure; imgray(path_map_c);
    %save([d_root 'predictions/prob_paths/' zerostr(ii,2) '_path.mat'], 'group_map');    
    
end
%%
d_root = [asymmetryroot 'data/retinograms/dRIVE/test/'];
mkdir([d_root 'predictions/prob_paths/']);
for ii = 1:20
    %load ret and take green channel
    ret = load_uint8([d_root 'images_extended/' zerostr(ii,2) '_test_ext']);
    ret = ret(:,:,2);

    f_mask = u_load([d_root 'foveal_masks/' zerostr(ii,2) '_test_f_mask.mat']);
    
    %compute 1st deriv.
    [g1d_r g1d_o] = gaussian_1st_derivative_gradient2(ret, 1.5);

    %compute 1st deriv.
    [g2d_r g2d_o] = gaussian_2nd_derivative_line(ret, 1.5);

    %Compute ridges
    ridge_map = mb_curvature_sign_change(g1d_r, g2d_r, g2d_o, 0, 0) > 0;
    ridge_map(~f_mask) = 0;
    
    %Compute convex sets
    group_map = staal_ridge_group(ridge_map, g2d_o);
    
    %Compute convex set regions
    [py px] = find(f_mask & ~group_map);
    [ry rx] = find(group_map > 0);
    
    set_nums = griddata(rx, ry, group_map(group_map > 0), px, py, 'nearest');
    
    set_map = group_map;
    set_map(f_mask & ~group_map) = set_nums;
    figure; imgray(set_map); colormap lines;
    
    %save outputs
    save([d_root 'predictions\ridges\g2d\' zerostr(ii,2) '_ridge.mat'], 'ridge_map', 'group_map', 'set_map');
    
    
    %figure; imgray(ridge);
end
%%
d_root = [asymmetryroot 'data/retinograms/dRIVE/test/'];
for ii = 1:20
    load([d_root 'predictions/prob_paths/' zerostr(ii,2) '_path.mat'], 'path_map_c');
    v_mask = u_load([d_root 'vessel_masks/' zerostr(ii,2) '_test_v_mask.mat']);
    figure; 
    subplot(2,2,1); imgray(path_map_c);
    subplot(2,2,2); imgray(v_mask);
    subplot(2,2,1); imgray(bwmorph(path_map_c, 'thin);
    subplot(2,2,2); imgray(v_mask);
end
%%
max_val = 0;
for ii = 1:20
    %path_map = u_load([d_root 'predictions/prob_paths/no_edge/' zerostr(ii,2) '_path_c.mat']);
    %path_map = u_load([d_root 'predictions/prob_paths/no_j/5e5/' zerostr(ii,2) '_path_map_c.mat']);
    path_map = u_load([d_root 'predictions/prob_paths/no_edge/' zerostr(ii,2) '_path.mat']);
    max_val = max(max_val, max(log(path_map(:))));
end
    
%%
d_root = [asymmetryroot 'data/retinograms/dRIVE/test/'];

tp_counts = zeros(20, 102);
fp_counts = zeros(20, 102);
t_counts = zeros(20, 1);
f_counts = zeros(20, 1);

centre_vessels = 0;
discard_edge = 0;
nms = 0;
data_type = 11;
for ii = 1:20
    
    switch data_type
        case 1
            label = 'Prob paths, centre-line, no junctions';
            load([d_root 'predictions/prob_paths/no_j/5e5/' zerostr(ii,2) '_path_map_c.mat'], 'path_map_c');
            vessel_prob = log(path_map_c) ./ max_val; %max(log(path_map_c(:)));%prctile(path_map_c(:), 99.5);
            vessel_prob(vessel_prob == -inf) = 0;
            
        case 2
            label = 'Raw classifications, centre-line';
            vessel_prob = load_uint8([d_root 'predictions/centre_detection/dt/rf_3/' zerostr(ii,2) '_test_ext_class.mat']);
            
        case 3
            label = 'Prob paths, centre-line, junctions';
            load([d_root 'predictions/prob_paths/j/' zerostr(ii,2) '_path_c.mat'], 'path_map_c');
            vessel_prob = log(path_map_c) ./ max(log(path_map_c(:)));%prctile(path_map_c(:), 99.5);
            vessel_prob(vessel_prob == -inf) = 0;
            
        case 4
            label = 'Raw classifications, all vessels';
            vessel_prob = load_uint8([d_root 'predictions/detection/dt/rf_3/' zerostr(ii,2) '_test_ext_class.mat']);
            
        case 5
            label = 'Prob paths, all vessels, no junctions';
            load([d_root 'predictions/prob_paths/no_j/' zerostr(ii,2) '_path.mat'], 'path_map');
            vessel_prob = log(path_map) ./ max(log(path_map(:)));%prctile(path_map_c(:), 99.5);
            vessel_prob(vessel_prob == -inf) = 0;
            
        case 6
            label = 'Prob paths, all vessels, junctions';
            load([d_root 'predictions/prob_paths/j/' zerostr(ii,2) '_path.mat'], 'path_map');
            vessel_prob = log(path_map) ./ max(log(path_map(:)));%prctile(path_map_c(:), 99.5);
            vessel_prob(vessel_prob == -inf) = 0;
        case 7
            label = 'Prob paths, centre-line, no junctions';
            load([d_root 'predictions/prob_paths/no_j/10e5/' zerostr(ii,2) '_path_map_c.mat'], 'path_map_c');
            vessel_prob = log(path_map_c) ./ max(log(path_map_c(:)));%prctile(path_map_c(:), 99.5);
            vessel_prob(vessel_prob == -inf) = 0;
        case 8
            label = 'Prob paths, centre-line, no junctions';
            load([d_root 'predictions/prob_paths/no_j/7.5e5/' zerostr(ii,2) '_path_map_c.mat'], 'path_map_c');
            vessel_prob = log(path_map_c) ./ max_val;%max(log(path_map_c(:)));%prctile(path_map_c(:), 99.5);
            vessel_prob(vessel_prob == -inf) = 0;
        case 9
            label = 'Prob paths, centre-line, no edge';
            load([d_root 'predictions/prob_paths/no_edge/' zerostr(ii,2) '_path_c.mat'], 'path_map_c');
            vessel_prob = log(path_map_c) ./ max_val;%max(log(path_map_c(:)));%prctile(path_map_c(:), 99.5);
            vessel_prob(vessel_prob == -inf) = 0;
        case 10
            label = 'Prob paths, all-line, no edge';
            load([d_root 'predictions/prob_paths/no_edge/' zerostr(ii,2) '_path.mat'], 'path_map');
            vessel_prob = log(path_map) ./ max(log(path_map(:)));%prctile(path_map_c(:), 99.5);
            vessel_prob(vessel_prob == -inf) = 0;
        case 11
            label = 'Raw classifications, all vessels';
            vessel_prob = load_uint8([d_root '/images/results/77469/' zerostr(ii,2) '_test_class.mat']);
        case 12
            label = 'Prob paths, all-line';
            load([d_root 'predictions/prob_paths/good/' zerostr(ii,2) '_path.mat'], 'path_map');
            vessel_prob = log(path_map) ./ max(log(path_map(:)));%prctile(path_map_c(:), 99.5);
            vessel_prob(vessel_prob == -inf) = 0;
        case 13
            label = 'Prob paths, all-line';
            vessel_prob = load_uint8([d_root '/images/results/77728/' zerostr(ii,2) '_test_class.mat']);
    end
    
    if nms
        [~, vessel_theta] = gaussian_2nd_derivative_line(vessel_prob, 1.5);
        vessel_prob = mb_non_maximal_supp(vessel_prob, vessel_theta);
    end
    
    v_mask = u_load([d_root 'vessel_masks\' zerostr(ii,2) '_test_v_mask.mat']);
    f_mask = u_load([d_root 'foveal_masks\' zerostr(ii,2) '_test_f_mask.mat']);
    if discard_edge
        d_mask = imerode(f_mask, strel('disk', 15));
        vessel_prob(f_mask & ~d_mask) = 0;
    end
    
    if centre_vessels
        v_mask = bwmorph(v_mask, 'skel', 'inf');
    end
    
    %Compute ROC counts for image
    [roc_pts auc tp_count fp_count] = calculate_roc_curve(vessel_prob(f_mask), v_mask(f_mask),(-1:100)/100);
    %[roc_pts auc tp_count fp_count] = calculate_roc_connected(vessel_prob, v_mask,(-1:100)/100, f_mask);
    
    %Increment total counts
    tp_counts(ii,:) = tp_count;
    fp_counts(ii,:) = fp_count;
    t_counts(ii) = sum(v_mask(f_mask));
    f_counts(ii) = sum(~v_mask(f_mask));

end

%Compute ROC points for complete set of data
roc_pts = [sum(fp_counts)' / sum(f_counts) sum(tp_counts)' / sum(t_counts)];

%Compute AUC for ROC curve
auc = sum( (roc_pts(2:end,1)-roc_pts(1:end-1,1)) .* roc_pts(1:end-1,2)) + ...
        0.5*sum( (roc_pts(2:end,1)-roc_pts(1:end-1,1)) .* (roc_pts(2:end,2)-roc_pts(1:end-1,2)) );
    
figure; 
plot(roc_pts(:,1), roc_pts(:,2), 'b'); axis([0 1 0 1]); axis equal; hold on;
plot(roc_pts(:,1), roc_pts(:,2), 'r.');
plot(roc_pts(51,1), roc_pts(51,2), 'gx');
title([ label ' AUC: ' num2str(auc)]);
%%
d_root = [asymmetryroot 'data/retinograms/dRIVE/test/'];

for ii = 1:20
    path_map_c = u_load([d_root 'predictions/prob_paths/no_j/' zerostr(ii,2) '_path_map_c.mat']);
    %path_map_cj = u_load([d_root 'predictions/prob_paths/j/' zerostr(ii,2) '_path_c.mat']);
    vessel_prob_c = load_uint8([d_root 'predictions/centre_detection/dt/rf_3/' zerostr(ii,2) '_test_ext_class.mat']);
    
    figure;
    a1 = subplot(1,2,1); imgray(log(path_map_c));
    a2 = subplot(1,2,2); imgray(log(path_map_c) ./ vessel_prob_c);
    linkaxes([a1 a2]);
end
%%
for ii = 1:20
    path_map_j = u_load([d_root 'predictions/prob_paths/j/' zerostr(ii,2) '_path.mat']);
    vessel_prob = load_uint8([d_root 'predictions/detection/dt/rf_3/' zerostr(ii,2) '_test_ext_class.mat']);
    f_mask = u_load([d_root 'foveal_masks\' zerostr(ii,2) '_test_f_mask.mat']);
    vessel_prob(~f_mask) = 0;
    
    figure;
    a1 = subplot(1,2,1); imgray(log(path_map_j));
    a2 = subplot(1,2,2); imgray(vessel_prob);
    linkaxes([a1 a2]);
end
%%
d_root = [asymmetryroot 'data/retinograms/dRIVE/test/'];
for ii = 1:20
%     path_map_c = u_load([d_root 'predictions/prob_paths/no_j/' zerostr(ii,2) '_path_map_c.mat']);
%     path_map_cj = u_load([d_root 'predictions/prob_paths/j/' zerostr(ii,2) '_path_c.mat']);
    
    path_map = u_load([d_root 'predictions/prob_paths/no_j/' zerostr(ii,2) '_path.mat']);
%     path_map_j = u_load([d_root 'predictions/prob_paths/j/' zerostr(ii,2) '_path.mat']);
    
%     log_path_c = log(path_map_c) / max(log(path_map_c(:))); log_path_c(log_path_c == -inf) = 0;
%     log_path_cj = log(path_map_cj) / max(log(path_map_cj(:))); log_path_cj(log_path_cj == -inf) = 0;
    log_path = log(path_map) / max(log(path_map(:))); log_path(log_path == -inf) = 0;
%     log_path_j = log(path_map_j) / max(log(path_map_j(:))); log_path_j(log_path_j == -inf) = 0;
    
%     vessel_prob = load_uint8([d_root 'predictions/detection/dt/rf_3/' zerostr(ii,2) '_test_ext_class.mat']);
%     vessel_prob_c = load_uint8([d_root 'predictions/centre_detection/dt/rf_3/' zerostr(ii,2) '_test_ext_class.mat']);
%     
%     f_mask = u_load([d_root 'foveal_masks\' zerostr(ii,2) '_test_f_mask.mat']);
%     vessel_prob(~f_mask) = 0;
%     vessel_prob_c(~f_mask) = 0;
    
    write_im_from_colormap(path_map, [d_root 'temp_figs/prob_paths/all/' zerostr(ii,2) '_prob_paths.png'], gray(256));
    write_im_from_colormap(log_path, [d_root 'temp_figs/log_paths/all/' zerostr(ii,2) '_prob_paths.png'], gray(256));
    %write_im_from_colormap(path_map_c, [d_root 'temp_figs/prob_paths/centre/' zerostr(ii,2) '_prob_paths.png'], gray(256));
%     write_im_from_colormap(log_path_c, [d_root 'temp_figs/log_paths/centre/' zerostr(ii,2) '_prob_paths.png'], gray(256));
    
    %write_im_from_colormap(path_map_j, [d_root 'temp_figs/prob_paths/all/j/' zerostr(ii,2) '_prob_paths.png'], gray(256));
%     write_im_from_colormap(log_path_j, [d_root 'temp_figs/log_paths/all/j/' zerostr(ii,2) '_prob_paths.png'], gray(256));
    %write_im_from_colormap(path_map_cj, [d_root 'temp_figs/prob_paths/centre/j/' zerostr(ii,2) '_prob_paths.png'], gray(256));
%     write_im_from_colormap(log_path_cj, [d_root 'temp_figs/log_paths/centre/j/' zerostr(ii,2) '_prob_paths.png'], gray(256));
    
    %write_im_from_colormap(vessel_prob, [d_root 'temp_figs/raw_probs/all/' zerostr(ii,2) '_raw_probs.png'], gray(256));
    %write_im_from_colormap(vessel_prob_c, [d_root 'temp_figs/raw_probs/centre/' zerostr(ii,2) '_raw_probs.png'], gray(256));
end
%%
d_root = [asymmetryroot 'data/retinograms/dRIVE/test/'];

for ii = 1:20
    path_map_c = u_load([d_root 'predictions/prob_paths/no_j/' zerostr(ii,2) '_path_map_c.mat']);
    path_map_e = u_load([d_root 'predictions/prob_paths/no_edge/' zerostr(ii,2) '_path_c.mat']);
    
    figure;
    a1 = subplot(1,2,1); imgray(log(path_map_c));
    a2 = subplot(1,2,2); imgray(log(path_map_e));
    linkaxes([a1 a2]);
end
%%
NUM_TREES=10 PREDICTION_TYPE="orientation" IMAGE_DIR="retinograms/DRIVE/training" ORI_DIR="orientations" SAMPLING_METHOD="generate_vessel_data" WIN_SIZE=3 DO_UBOUND=0 SPLIT_CRITERION="ssq" VAR_CRITERION="ssq" END_CUT_MIN=0 DECOMP_TYPE="g2d" SIGMA_RANGE="[0.5 1 2 4 8]" qsub -t 1-20 -V -l twoday matlab_code/trunk/hydra/cuc/build_vessel_predictor.sh
FOREST_JOB="'48947'" FOREST_DIR="line_orientation_rfs" qsub -V -l short matlab_code/trunk/hydra/cuc/combine_hydra_rfs.sh
NUM_JOBS=20 FOREST_JOB="'48947'" TEST_IMAGE_DIR="'retinograms/DRIVE/test/images_extended'" MASK_DIR="retinograms/DRIVE/test/foveal_masks" FOREST_DIR="line_orientation_rfs" qsub -l short -t 1 -V matlab_code/trunk/hydra/cuc/classify_image_set.sh
qstat

