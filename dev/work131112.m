apex_gt_dir = 'C:\isbe\nailfold\data\rsa_study\test\apex_gt\';
test_list = dir([apex_gt_dir '*.mat']);
num_images = length(test_list);
%%

for i_im = 1:20
    
    im_num = test_list(i_im).name(1:6);
    
    load([apex_gt_dir im_num '_gt.mat']);

    all_xy = apex_xy;

    num_candidates = size(all_xy,1);
    
    if num_candidates < 10
        continue;
    end

    [distal_idx] = select_distal_candidates(all_xy, 10, 1);

    valid_apexes = is_distal;% & (num_apex_markers == num_im_markers);
    
    plot(apex_xy(valid_apexes,1), apex_xy(valid_apexes,2), 'g+');
    %plot(apex_xy(is_undefined,1), apex_xy(is_undefined,2), 'bs');
 
end
%%
for i_im = 61:90
    im_name = im_list(i_im).name(1:6);
    
    load([apex_gt_dir im_name '_gt']);    
    
    if num_im_markers == 1
        continue;
    end
    
    load([vessel_centres_dir im_name '_vc'],...
        'vessel_centre');
    load([apex_map_dir 'local_maxima\' im_name '_candidates'],...
        'candidate_xy', 'candidate_scores');
    
    valid_apexes = is_distal & (num_apex_markers > 1);
    
    figure;
    plot(vessel_centre.x, vessel_centre.y, 'b.', 'markersize', 2);
    axis image ij; hold all;
    title(im_name);
    
    plot(apex_xy(valid_apexes,1), apex_xy(valid_apexes,2), 'g', 'linewidth', 2);
    plot(apex_xy(is_distal,1), apex_xy(is_distal,2), 'gx');
    
    num_candidates = size(candidate_xy,1);
    
    confirmed_xy = [];
    confirmed_scores = [];
    potential_xy = [];
    potential_scores = [];
    secondary_xy = [];
    secondary_scores = [];
    
    
    valid_candidates = candidate_scores > 0.15;
    
    
    if any(valid_candidates)
%         [distal_idx] = select_distal_candidates(candidate_xy(valid_candidates,:),...
%             candidate_scores(valid_candidates), 10, [], 0);
        
        [distal_idx] = select_distal_candidates(candidate_xy(valid_candidates,:),...
            [], 10, [], 0);
    
        potential_xy = candidate_xy(valid_candidates,:);
        potential_xy = potential_xy(distal_idx,:);
        potential_xy = sortrows(potential_xy);
        
        plot(candidate_xy(valid_candidates,1), candidate_xy(valid_candidates,2), 'r*');
        plot(candidate_xy(~valid_candidates,1), candidate_xy(~valid_candidates,2), 'y*');
        
        plot(potential_xy(:,1), potential_xy(:,2), 'r-', 'linewidth', 2);
        
    end
    
    
%     if num_candidates > stage1_max
%         
%         max_candidates = min(num_candidates, 20);
%         [distal_idx] = select_distal_candidates(candidate_xy(1:max_candidates,:), 10, [], 0);
%         
%         confirmed_xy = candidate_xy(distal_idx,:);
%         confirmed_scores = candidate_scores(distal_idx,:);
%         
%         potential_xy = candidate_xy(stage2_max+1:end,:);
%         potential_scores = candidate_scores(stage2_max+1:end,:);
%         
%         
%         if num_candidates > stage2_max                       
%             [distal_idx] = evaluate_distal_candidates(confirmed_xy, potential_xy);   
%             
%             %plot(potential_xy(distal_idx,1), potential_xy(distal_idx,2), 'bx');
%             %plot(potential_xy(~distal_idx,1), potential_xy(~distal_idx,2), 'mv');
%             
%             secondary_xy = potential_xy(distal_idx,:);
%             secondary_scores = potential_scores(distal_idx,:);
%             
%             
%         end
%         
%         candidate_xy = [confirmed_xy; secondary_xy];
%         candidate_scores = [confirmed_scores; secondary_scores];
%         
%         [distal_idx] = select_distal_candidates(candidate_xy, 10, [], 0);
%         
%         candidate_xy = candidate_xy(distal_idx,:);
%         candidate_scores = candidate_scores(distal_idx,:);
%         
%         %plot(candidate_xy(:,1), candidate_xy(:,2), 'go');
%     end   
%     
%     save([apex_map_dir 'distal_maxima\' im_name '_maxima'],...
%         '*_xy', '*_scores');
end

%%
for i_im = 1:20%num_images;%[ 6 40 92 202 219 223 255 266 323 344 494 517 546]%
    im_num = test_list(i_im).name(1:6);
    load([apex_gt_dir im_num '_gt.mat'], 'distal_xy', 'distal_width', 'non_distal_xy', 'undefined_xy');
    load(['C:\isbe\nailfold\data\rsa_study\test\vessel_centres\' im_num '_vc.mat']);
    
    distal_row_x = apex_xy(is_distal,1) - ncols/2;
    if length(distal_row_x) < 5
        continue;
    end
    
    distal_row_y = apex_xy(is_distal,2);
    
    include_pts = vessel_centre_discards ~= 1;
    
    centre_vecs = angle(vessel_centre_ori(include_pts)) / 2;
    centre_vecs = [cos(centre_vecs) sin(centre_vecs)];

    vx = vessel_centre_x(include_pts) - ncols/2;    
    vy = vessel_centre_y(include_pts) - nrows/2;   
    %vw = vessel_centre_corr(include_pts) .* vessel_centre_prob(include_pts);
    %vw(vw < 0) = 0;
    vw = vessel_centre_prob(include_pts);
    
    
    pp = polyfit(distal_row_x, distal_row_y, 1);
    lobf_x = [-ncols ncols]/2;
    lobf_y = pp(1)*lobf_x + pp(2);  
    
    targets(i_im,:) = pp;
    
    do_plot_i = do_plot && i_im <= 10;
    if do_plot_i
        figure;
        a1 = subplot(3,4,1);
        axis image ij; axis([-ncols ncols -nrows nrows]/2); hold on; 
        plot(distal_xy(:,1) - ncols/2, distal_xy(:,2) - nrows/2, 'rx');
        plot(undefined_xy(:,1) - ncols/2, undefined_xy(:,2) - nrows/2, 'bx');
        plot(lobf_x, lobf_y, 'g');
    
        plot_num = 2;
    end
    
    max_angle_sum = 0;
    max_angle = 0;
    offset_mean = 0;
    for i_ra = 1:num_angles
        ref_angle = ref_angles(i_ra);

        projection_weights = abs(centre_vecs * [cos(ref_angle+pi/2); sin(ref_angle)+pi/2]);
        vyi = [vx vy]*[-sin(ref_angle); cos(ref_angle)];
        
        valid_pts = (vyi > -200) & (vyi < 200);
        
        [y_hist_angle_weighted] = compute_weighted_histogram(vyi(valid_pts), vw(valid_pts) .* projection_weights(valid_pts), bin_centres);       
        %feature_vectors(i_im, i_ra, :) = y_hist_angle_weighted(2:end-1) / sum(vessel_centre_prob);
        
        %angle_sum_i = sum(vw .* projection_weights);
        angle_sum_i = sum(y_hist_angle_weighted);
        
        display(angle_sum_i);
        offset_mean_i = sum(vyi(valid_pts) .* vw(valid_pts) .* projection_weights(valid_pts)) / angle_sum_i;
        if angle_sum_i > max_angle_sum
            max_angle_sum = angle_sum_i;
            max_angle = ref_angle;
            offset_mean = offset_mean_i;
        end
        
        max_angles(i_im) = max_angle;
        offset_means(i_im) = offset_mean;
        
        if do_plot_i
            subplot(3,4,plot_num);
            bar(bin_centres, y_hist_angle_weighted / sum(vw(valid_pts)));
            set(gca, 'ylim', [0 0.25]);
            plot_num = plot_num+1;
            title(['(' num2str(sum(y_hist_angle_weighted)) ', ' num2str(sum(y_hist_angle_weighted)) ')']);
        end
        title(a1, ['Max angle = ' num2str(round(180*max_angle/pi)), ' \theta = ' num2str(num2str(round(180*atan(pp(1))/pi)))]);
        xlabel(a1, ['Offset = ' num2str(round(offset_mean)) ' y = ' num2str(round(pp(2)))]);
    end
    
    %load(['C:\isbe\nailfold\data\rsa_study\test\images\' im_num '.mat']);
    %figure; imgray(nailfold);
end
%%
im_list = dir('C:\isbe\nailfold\data\rsa_study\test\images\*.mat');

apex_map_dir = 'C:\isbe\nailfold\data\rsa_study\test\apex_maps\frog\';
apex_gt_dir = [nailfoldroot 'data/rsa_study/test/apex_gt/'];
vessel_centres_dir = 'C:\isbe\nailfold\data\rsa_study\test\vessel_centres\';

%%
for i_im = 12:20%length(im_list);
    im_name = im_list(i_im).name(1:6);  
    
    load([apex_map_dir 'distal_maxima\' im_name '_maxima']);
    load([apex_gt_dir im_name '_gt']);    
      
    if num_im_markers == 1
        continue;
    end
    
    load([vessel_centres_dir im_name '_vc'],...
        'vessel_centre');
    
    valid_apexes = is_distal & (num_apex_markers > 1);
    
    figure;
    plot(vessel_centre.x, vessel_centre.y, 'b.', 'markersize', 2);
    axis image ij; hold all;
    
    plot(apex_xy(valid_apexes,1), apex_xy(valid_apexes,2), 'g', 'linewidth', 2);
    plot(apex_xy(is_distal,1), apex_xy(is_distal,2), 'gx');
    
    if ~isempty(confirmed_xy)
        confirmed_xy = sortrows(confirmed_xy);
        plot(confirmed_xy(:,1), confirmed_xy(:,2), 'r', 'linewidth', 2);
        plot(confirmed_xy(:,1), confirmed_xy(:,2), 'rx');
    end    
    
    load([apex_map_dir 'local_maxima\' im_name '_candidates'],...
        'candidate_xy', 'candidate_scores');
    
    plot(candidate_xy(:,1), candidate_xy(:,2), 'mv');
    
    if ~isempty(potential_xy)
        rejected_xy = setdiff(candidate_xy, potential_xy, 'rows');
        
        if ~isempty(rejected_xy)
            plot(rejected_xy(:,1), rejected_xy(:,2), 'y^');
            plot(rejected_xy(:,1), rejected_xy(:,2), 'yv');
        end
    end
        
    if ~isempty(secondary_xy)
        plot(secondary_xy(:,1), secondary_xy(:,2), 'm^');
    end
    
end
%%

figure; hold on;
for i_im = 1:length(im_list)
    im_name = im_list(i_im).name(1:6);
    
    load([apex_gt_dir im_name '_gt']);    
    
    if num_im_markers == 1
        continue;
    end
    
    load([apex_map_dir 'local_maxima\' im_name '_candidates'],...
        'candidate_xy', 'candidate_scores');
    
    valid_apexes = is_distal & (num_apex_markers > 1);
    
    valid_candidates = candidate_scores > 0.127;
    
    plot(sum(valid_apexes), sum(valid_candidates), 'bx');
end
%%
im_list = dir('C:\isbe\nailfold\data\rsa_study\test\images\*.mat');
fov_mask_dir = 'C:\isbe\nailfold\data\rsa_study\test\fov_masks\';
im_rot_dir = 'C:\isbe\nailfold\data\rsa_study\test\image_rotations\';

ori_bins = -45:45; %linspace(0, pi, 180);

for i_im = 1:length(im_list);
    im_name = im_list(i_im).name(1:6);
    
    f_mask = u_load([fov_mask_dir im_name '_f_mask.mat']);
    
    f_mask4 = imresize(f_mask, 1/8, 'nearest');
    
    bad_mask = imopen(f_mask4, strel('disk', 55));
    
    bad = any(bad_mask(:));
     
    f_border = f_mask & ~imerode(f_mask, strel('disk', 2));
    
    [~, border_ori_map] = gaussian_1st_derivative_gradient(f_mask, 4);
    
    border_ori = border_ori_map(f_border);
    border_ori(border_ori < -pi/4) = [];
    border_ori(border_ori > pi/4) = [];
    border_ori = 180*border_ori/pi;
    
    counts = hist(border_ori, ori_bins);
    [~, max_idx] = max(counts);
    max_ori = ori_bins(max_idx);
    
    r_mask = imrotate(f_mask, -max_ori, 'nearest', 'loose');
    [nrows ncols] = size(f_mask);
    [nrowsr ncolsr] = size(r_mask);

    depth = zeros(ncolsr, 1);
    centres_r = zeros(ncolsr, 2);
    centres_r(:,1) = 1:ncolsr;
    for i_col = 1:ncolsr
        depth(i_col) = sum(r_mask(:,i_col));
        centres_r(i_col,2) = mean(find(r_mask(:,i_col)));
    end

    start_col = find(depth > 470, 1);
    end_col = find(depth > 470, 1, 'last');
    
    centres_r = centres_r(start_col:end_col,:);
    
    [pp, s, mu] = polyfit(centres_r(:,1), centres_r(:,2), 5);
    centres_r(:,2) = polyval(pp, centres_r(:,1), s, mu);
    
    rot_mat = [cosd(max_ori) -sind(max_ori); sind(max_ori) cosd(max_ori)];

    centres = bsxfun(@plus,  bsxfun(@minus, centres_r, [ncolsr/2 nrowsr/2])*rot_mat, [ncols nrows]/2);
    
    save([im_rot_dir im_name '.mat'], 'centres_r', 'rot_mat', 'ncolsr', 'nrowsr', 'ncols', 'nrows', 'bad');
    
end
%%
im_list = dir('C:\isbe\nailfold\data\rsa_study\test\images\*.mat');
fov_mask_dir = 'C:\isbe\nailfold\data\rsa_study\test\fov_masks\';
apex_gt_dir = 'C:\isbe\nailfold\data\rsa_study\test\apex_gt\';
aam_dir = 'C:\isbe\nailfold\data\rsa_study\test\aam\';
im_rot_dir = 'C:\isbe\nailfold\data\rsa_study\test\image_rotations\';
candidates_dir = 'C:\isbe\nailfold\data\rsa_study\test\apex_maps\frog\post_merged\';
ori_bins = -45:45; %linspace(0, pi, 180);

displacement_hist_bins = -400:10:400;
max_candidates = 5:5:25;%2:2:20;

apex_displacement_hist = zeros(length(max_candidates), length(displacement_hist_bins));

for i_im = 1:length(im_list);
    im_name = im_list(i_im).name(1:6);
    
    load([apex_gt_dir im_name '_gt']);    
    load([im_rot_dir im_name '.mat'], 'centres_r', 'rot_mat', 'ncolsr', 'nrowsr', 'ncols', 'nrows', 'bad');
    
    if num_im_markers == 1 || bad
        continue;
    end
    
    
    load([candidates_dir im_name '_candidates'],...
        'candidate_xy', 'candidate_scores');
    [candidate_scores cs_i] = sort(candidate_scores, 'descend');
    candidate_xy = candidate_xy(cs_i,:);
    
    if length(candidate_xy) < 3
        continue;
    end
    
    valid_apexes = is_distal & (num_apex_markers > 1);
    valid_candidates = candidate_scores > 0.3;
    candidate_xy = candidate_xy(valid_candidates,:);
    
    apex_xyr = bsxfun(@plus, bsxfun(@minus, apex_xy, [ncols nrows]/2)*rot_mat', [ncolsr nrowsr]/2);
    candidate_xyr = bsxfun(@plus, bsxfun(@minus, candidate_xy, [ncols nrows]/2)*rot_mat', [ncolsr nrowsr]/2);
    
    [pp, s, mu] = polyfit(centres_r(:,1), centres_r(:,2), 5);
    
    
    
    apex_polyfit = polyval(pp, apex_xyr(:,1), s, mu);
    apex_displacements = apex_xyr(:,2) - apex_polyfit;
    
        
    candidate_polyfit = polyval(pp, candidate_xyr(:,1), s, mu);
    candidate_displacements = candidate_xyr(:,2) - candidate_polyfit;
    
    for ii = 1:length(max_candidates)
        num_to_select = min(size(candidate_xy,1), max_candidates(ii));
        
        apex_displacement_hist(ii,:) = apex_displacement_hist(ii,:) +...
            hist(apex_displacements(valid_apexes) - median(candidate_displacements(1:num_to_select)), displacement_hist_bins); 
    end

    if i_im <= 0
        candidate_xyrp = sortrows(candidate_xyr(1:num_to_select,:));
        
        figure; axis equal ij; hold on;
        plot(centres_r(:,1), centres_r(:,2));
        plot(apex_xyr(valid_apexes,1), apex_xyr(valid_apexes,2), 'g', 'linewidth', 2);
        plot(candidate_xyrp(:,1), candidate_xyrp(:,2), 'r', 'linewidth', 2);
        plot(candidate_xyrp(:,1), candidate_xyrp(:,2), 'bx', 'linewidth', 2);
        plot(centres_r(:,1), centres_r(:,2) + mean(candidate_displacements), 'c');
        
        plot(candidate_xyr(1:num_to_select,1), candidate_polyfit, 'gx');
        plot(apex_xyr(:,1), apex_polyfit, 'rx');
    end
end

figure; hold all;
for ii = 1:length(max_candidates)
    plot(displacement_hist_bins, cumsum(apex_displacement_hist(ii,:)) / sum(apex_displacement_hist(ii,:)));
end
legend;
plot(displacement_hist_bins([1 end]), [0.95 0.95], 'k');
plot(displacement_hist_bins([1 end]), [0.05 0.05], 'k');
%%

max_candidates = 15;

total_hits = 0;
total_misses = 0;
total_kept = 0;
total_discarded = 0;

fixed_distal_detections = 0;
fixed_non_distal_detections = 0;
fixed_undefined_detections = 0;  
false_positives = 0;

load('C:\isbe\nailfold\data\rsa_study\test\results\post_merged_smooth_20131122T142242.mat');

figure; hold all;
xlabel('Displacement to mid-line');
ylabel('RF offset score');

displacement_hist_bins = -400:20:400;
rf_offset_bins = linspace(0, 5, 50);
bins = {displacement_hist_bins, rf_offset_bins};

fp_dist = zeros(length(displacement_hist_bins), length(rf_offset_bins));
tp_dist = zeros(length(displacement_hist_bins), length(rf_offset_bins));

all_cans = zeros(0,4);
all_labels = zeros(0,1);

for i_im = 1:length(im_list);
    im_name = im_list(i_im).name(1:6);
    
    load([apex_gt_dir im_name '_gt']);    
    load([im_rot_dir im_name '.mat'], 'centres_r', 'rot_mat', 'ncolsr', 'nrowsr', 'ncols', 'nrows', 'bad');
    
    if bad %num_im_markers == 1 || bad
        continue;
    end
    
    
    load([candidates_dir im_name '_candidates'],...
        'candidate_xy', 'candidate_scores');
    
    [candidate_scores cs_i] = sort(candidate_scores, 'descend');
    candidate_xy = candidate_xy(cs_i,:);

    if length(candidate_xy) < 3
        continue;
    end
    
    fid = fopen([aam_dir im_name '\model_qualities.txt']);
    q_txt = textscan(fid, '%s %f', 'delimiter', ':');
    fclose(fid);
    model_q = q_txt{2}; clear q_txt;
    
    valid_apexes = is_distal & (num_apex_markers > 1);
    valid_candidates = candidate_scores > 0.3;
    
    apex_xyr = bsxfun(@plus, bsxfun(@minus, apex_xy, [ncols nrows]/2)*rot_mat', [ncolsr nrowsr]/2);
    candidate_xyr = bsxfun(@plus, bsxfun(@minus, candidate_xy, [ncols nrows]/2)*rot_mat', [ncolsr nrowsr]/2);
    
    [pp, s, mu] = polyfit(centres_r(:,1), centres_r(:,2), 5); 
    
    apex_polyfit = polyval(pp, apex_xyr(valid_apexes,1), s, mu);
    apex_displacements = apex_xyr(valid_apexes,2) - apex_polyfit;
    
    candidate_polyfit = polyval(pp, candidate_xyr(:,1), s, mu);
    candidate_displacements = candidate_xyr(:,2) - candidate_polyfit;
    
    num_to_select = min(sum(candidate_scores > 0.3,1), max_candidates);
    offset = mean(candidate_displacements(1:num_to_select)); 
    total_displacement = mean(abs(candidate_displacements(1:num_to_select)));
    
    upper_lim = offset - 135; 
    lower_lim = offset + 80;
    
    hits = (apex_displacements > upper_lim) & (apex_displacements < lower_lim);
    total_hits = total_hits + sum(hits);
    total_misses = total_misses + sum(~hits);
    
    kept = (candidate_scores > 0.3) & ...
        ((candidate_displacements > upper_lim) | (candidate_scores > 0.8) ) &...
        (candidate_displacements < lower_lim);
    
    potential = ~kept & (candidate_scores > 0.8);
    [reclaimed] = evaluate_distal_candidates(candidate_xy(kept,:), ...
        candidate_xy(potential,:));
    
    kept(potential) = reclaimed;
    
    total_kept = total_kept + sum(kept);
    total_discarded = total_discarded + sum(~kept);    
    
    detected_idx = detections{i_im,3}(kept);
    true_detected_idx = ...
        unique(detected_idx(detected_idx > 0));

    included_distal_detections = ...
        valid_apexes(true_detected_idx)  & ...
        is_distal(true_detected_idx);

    included_non_distal_detections = ...
        valid_apexes(true_detected_idx)  & ...
        is_non_distal(true_detected_idx);

    included_undefined_detections = ...
        valid_apexes(true_detected_idx)  & ...
        is_undefined(true_detected_idx);

    fixed_distal_detections = fixed_distal_detections + ...
        sum(included_distal_detections);

    fixed_non_distal_detections = fixed_non_distal_detections + ...
        sum(included_non_distal_detections);

    fixed_undefined_detections = fixed_undefined_detections + ...
        sum(included_undefined_detections);

    false_positives = false_positives + sum(~detected_idx);
    
    tp = detections{i_im,3} > 0;
    tp(tp) = is_distal(detections{i_im,3}(tp));    
    fp = ~detections{i_im,3};
    nd = detections{i_im,3} > 0;
    nd(nd) = is_non_distal(detections{i_im,3}(nd));    
    %fp(~fp) = is_non_distal(detections{i_im,3}(~fp));
    
    tp_dist = tp_dist + hist3([candidate_displacements(tp) - offset candidate_scores(tp)],...
        bins); 
    fp_dist = fp_dist + hist3([candidate_displacements(fp) - offset candidate_scores(fp)],...
        bins); 
    
    all_cans = [all_cans;...
        candidate_displacements(tp) - offset...
        candidate_scores(tp)...
        (candidate_displacements(tp) - offset) / total_displacement ...
        model_q(tp)]; %#ok
    all_cans = [all_cans;...
        candidate_displacements(fp) - offset...
        candidate_scores(fp)...
        (candidate_displacements(fp) - offset) / total_displacement ...
        model_q(fp)]; %#ok
    all_cans = [all_cans;...
        candidate_displacements(nd) - offset...
        candidate_scores(nd)...
        (candidate_displacements(nd) - offset) / total_displacement ...
        model_q(nd)]; %#ok  
    
    all_labels = [all_labels; zeros(sum(tp),1); ones(sum(fp),1);  2*ones(sum(nd),1);]; %#ok
    
    
    plot(candidate_displacements(fp) - offset,...
        candidate_scores(fp), 'rx');
    plot(candidate_displacements(tp) - offset,...
        candidate_scores(tp), 'g+');
    plot(candidate_displacements(nd) - offset,...
        candidate_scores(nd), 'bo');
    
    if i_im <= 0
        figure; 
        subplot(2,1,1); axis equal ij; hold on;
        title(im_name);
        plot(centres_r(:,1), centres_r(:,2) + offset, 'c');
        plot(centres_r(:,1), centres_r(:,2) + offset - 130, 'c--');
        plot(centres_r(:,1), centres_r(:,2) + offset + 100, 'c--');

        plot(candidate_xyr(:,1), candidate_xyr(:,2), 'bx');
        plot(candidate_xyr(1:num_to_select,1), candidate_xyr(1:num_to_select,2), 'gs');
        plot(candidate_xyr(candidate_scores > 0.8,1), candidate_xyr(candidate_scores > 0.8,2), 'y*');
        plot(apex_xyr(valid_apexes,1), apex_xyr(valid_apexes,2), 'ro');
        
        subplot(2,1,2); axis equal ij; hold on;
        title(im_name);
        plot(centres_r(:,1), centres_r(:,2) + offset, 'c');
        plot(centres_r(:,1), centres_r(:,2) + offset - 130, 'c--');
        plot(centres_r(:,1), centres_r(:,2) + offset + 100, 'c--');

        plot(candidate_xyr(kept,1), candidate_xyr(kept,2), 'bx');
        plot(apex_xyr(is_distal,1), apex_xyr(is_distal,2), 'go');
        plot(apex_xyr(valid_apexes,1), apex_xyr(valid_apexes,2), 'ro');
    end

end
%%
keep = (all_cans(:,2) > 0.3) & ...
       ((all_cans(:,1) > -140) | (all_cans(:,2) > 0.8) ) &...
       (all_cans(:,1) < 90);
    
figure; hold on;

plot3(all_cans(all_labels & keep,1), all_cans(all_labels & keep,2), all_cans(all_labels & keep,4), 'gx');
plot3(all_cans(~all_labels & keep,1), all_cans(~all_labels & keep,2), all_cans(~all_labels & keep,4), 'r+');
xlabel('Displacement to image mid-line');
ylabel('RF offset score');
zlabel('AAM model fit');

%%
x = all_cans;
y = all_labels + 1;
discard = isnan(x(:,1));
x(discard,:) = [];
y(discard) = [];

x_pos = x(x(:,1) > 0,:);
y_pos = y(x(:,1) > 0,:);

x_neg = x(x(:,1) <= 0,:);
y_neg = y(x(:,1) <= 0,:);

figure; hold on;
plot(x_pos(y_pos==2,1), x_pos(y_pos==2,2), 'b+');
plot(x_pos(y_pos==1,1), x_pos(y_pos==1,2), 'rx');

[c,err,post,logl,str] = classify(x_pos,x_pos,y_pos,'quadratic');

K = str(1,2).const;
L = str(1,2).linear;
Q = str(1,2).quadratic;
% Plot the curve K + [x,y]*L + [x,y]*Q*[x,y]' = 0:
f = @(x,y) K + L(1)*x + L(2)*y ...
            + Q(1,1)*x.^2 + (Q(1,2)+Q(2,1))*x.*y + Q(2,2)*y.^2;
ezplot(f,[-400 400 0 5]);

figure; hold on;
plot(x_neg(y_neg==2,1), x_neg(y_neg==2,2), 'b+');
plot(x_neg(y_neg==1,1), x_neg(y_neg==1,2), 'rx');
[c,err,post,logl,str] = classify(x_neg,x_neg,y_neg,'quadratic');

K = str(1,2).const;
L = str(1,2).linear;
Q = str(1,2).quadratic;
% Plot the curve K + [x,y]*L + [x,y]*Q*[x,y]' = 0:
f = @(x,y) K + L(1)*x + L(2)*y ...
            + Q(1,1)*x.^2 + (Q(1,2)+Q(2,1))*x.*y + Q(2,2)*y.^2;
ezplot(f,[-400 400 0 5]);

%%
x = candidate_xyr(kept,1);
y = candidate_xyr(kept,2);
w = candidate_scores(kept);
w = w / mean(w);

[pp_o] = polyfit(x, y, 1);

modelFun = @(b,x) polyval(pp, s, mu);
yw = sqrt(w).*y;
modelFunw = @(b,x) sqrt(w).*modelFun(b,x);

[pp_w,rw,Jw,Sigmaw,msew] = nlinfit(x,yw,modelFunw,pp_o);

centres_w = centres_r;
[centres_w(:,2), deltaw] = nlpredci(modelFun, centres_w(:,1), pp_w, rw, 'cov', Sigmaw);
%%
tgt_x = candidate_xyr(kept,1);
tgt_y = candidate_xyr(kept,2); 

src_x = candidate_xyr(kept,1);
src_y = candidate_polyfit(kept);

[D, Z, T] = procrustes([tgt_x tgt_y], [src_x src_y], 'scaling', false, 'reflection', false);

figure; axis equal ij; hold on;
plot(tgt_x, tgt_y, 'rx');
plot(src_x, src_y, 'bx');
plot(Z(:,1), Z(:,2), 'gx');
plot(centres_w(:,1), centres_w(:,2), 'b');

centres_w = bsxfun(@plus, T.b * centres_r * T.T, T.c(1,:));


%%
for i_im = 1:30%length(im_list);
    
    im_name = im_list(i_im).name(1:6);
   
    f_mask = u_load([fov_mask_dir im_name '_f_mask.mat']);
    
    f_mask4 = imresize(f_mask, 1/8, 'nearest');
    
    bad_mask = imopen(f_mask4, strel('disk', 55));
    
    bad = any(bad_mask(:));
    
    figure;
    subplot(2,1,1); imgray(f_mask);
    if bad
        title('Image is naughty');
    else
        title('Image is ok!');
    end
    subplot(2,1,2); imgray(bad_mask);
end
%%
im_list = dir('C:\isbe\nailfold\data\rsa_study\test\images\*.mat');

apex_map_dir = 'C:\isbe\nailfold\data\rsa_study\test\apex_maps\frog\';
apex_gt_dir = [nailfoldroot 'data/rsa_study/test/apex_gt/'];
vessel_centres_dir = 'C:\isbe\nailfold\data\rsa_study\test\vessel_centres\';

%
for i_im = 1:30%length(im_list);
    im_name = im_list(i_im).name(1:6);  
    
    load([apex_map_dir 'local_maxima2\' im_name '_candidates.mat']);
    %load([apex_gt_dir im_name '_gt']);
    
    if length(candidate_scores) < 10
        continue;
    end
      
    [idx c sumd] = kmeans(candidate_xy, 6, 'EmptyAction', 'drop', 'Replicates', 100);
   
    figure; 
    subplot(2,1,1); axis equal ij; hold all;
    for i_k = 1:3; 
        plot(candidate_xy(idx==i_k,1), candidate_xy(idx==i_k,2), 'x'); 
    end
    title(im_name);
    
    subplot(2,1,2); axis equal ij; hold all;
    plot(c(:,1), c(:,2), 'rx');
    
end
%%
prob_dir = 'C:\isbe\nailfold\data\rsa_study\test\predictions\detection\rf_classification\257273\';
dist_thresh = 100;
patch_sz2 = 51;
connect_thresh = 0.5;
n_connect_pts = 20;

i_im = 3;
im_name = im_list(i_im).name(1:6);
    
load([candidates_dir im_name '_candidates'], 'candidate_xy', 'candidate_scores');
    
vessel_prob = u_load([prob_dir im_name '_pred.mat']);

[candidate_xy2, candidate_scores2] = post_merge_candidates(candidate_xy, candidate_scores, ...
    vessel_prob, dist_thresh, patch_sz2, connect_thresh, n_connect_pts, 1);
%%
%%
im_list = dir('C:\isbe\nailfold\data\rsa_study\test\images\*.mat');
fov_mask_dir = 'C:\isbe\nailfold\data\rsa_study\test\fov_masks\';
image_dir = 'C:\isbe\nailfold\data\rsa_study\test\images\';
im_rot_dir = 'C:\isbe\nailfold\data\rsa_study\test\image_rotations\';
candidates_dir = 'C:\isbe\nailfold\data\rsa_study\test\apex_maps\frog\post_merged\';

max_candidates = 15;

for i_im = 1:20%length(im_list);
    im_name = im_list(i_im).name(1:6);
    
    load([image_dir im_name '.mat']);
    load([fov_mask_dir im_name '_f_mask.mat']);
    load([candidates_dir im_name '_candidates'],...
        'candidate_xy', 'candidate_scores');   
    rot_struc = load([im_rot_dir im_name '.mat'], 'centres_r', 'rot_mat', 'ncolsr', 'nrowsr', 'ncols', 'nrows', 'bad');
    
    %[candidate_scores c_idx] = sort(candidate_scores);
    %candidate_xy = candidate_xy(c_idx,:);
    
    figure; imgray(nailfold);
    caxis([min(nailfold(fov_mask))+5 max(nailfold(fov_mask))-5]);
    
    if length(candidate_xy) < 3
        title([im_name ': very few candidates found, image may be of ungradeable quality']);
        continue;
    end
    
    if rot_struc.bad
        title({[im_name ': mosaic shape suggests frames were not properly registered.'],;...
            'No attempt made to find distal row. Vessel detections at 95% specificty shown'});
        kept = candidate_scores > 0.6;
        non_distal = false(size(kept));
        fell_at_the_last = [];
    else
        [kept non_distal intermediate_selections] = ...
            select_vessels_from_candidates(candidate_xy, candidate_scores, rot_struc);
        fell_at_the_last = intermediate_selections(:,3) & ~kept;
        
        title([im_name ': ' num2str(sum(kept)) ' distal apices found']);
    end

    if any(kept)
        plot(candidate_xy(kept,1), candidate_xy(kept,2), 'rx', 'markersize', 8);
    end
    if any(fell_at_the_last)
        plot(candidate_xy(fell_at_the_last,1), candidate_xy(fell_at_the_last,2), 'bo', 'markersize', 8);
    end
    if any(non_distal)
        plot(candidate_xy(non_distal,1), candidate_xy(non_distal,2), 'gx', 'markersize', 4);
    end
    if any(~kept & ~non_distal)
        plot(candidate_xy(~kept & ~non_distal,1), candidate_xy(~kept & ~non_distal,2), 'cx', 'markersize', 2);
    end

end

%%
select_vessels_from_candidates_batch('start_i', 1, 'end_i', 20);
%%

im_list = dir('C:\isbe\nailfold\data\rsa_study\test\images\*.mat');
prob_dir = 'C:\isbe\nailfold\data\rsa_study\test\predictions/detection/rf_classification/257273/';
candidates_dir = 'C:\isbe\nailfold\data\rsa_study\test\apex_maps\frog\post_merged\';
gt_dir = 'C:\isbe\nailfold\data\rsa_study\test\apex_gt\';

load('C:\isbe\nailfold\data\rsa_study\test\results\rf_offsets_20131118T120345.mat');

% i_im = 11;
% im_name = im_list(i_im).name(1:6); %'11839c';
load([candidates_dir im_name '_candidates'], 'candidate_xy', 'candidate_scores');
load([gt_dir im_name '_gt']);
%%    
for i_can = 1:20
    aam_can = u_load(['C:\isbe\nailfold\data\rsa_study\test\aam\' im_name '\apex' zerostr(i_can,4) '_aam.mat']);
    
    if i_can == 1
        nailfold = u_load(aam_can.image_path);       
        vessel_prob = u_load([prob_dir im_name '_pred.mat']);
    end
        
    if kept(i_can)
        pred_status = 'Kept :' ;
    else
        pred_status = 'Rejected :' ;
    end
    
    if detections{i_im,3}(i_can)       
        if is_distal(detections{i_im,3}(i_can))
            true_status = 'Distal: ';
        elseif is_non_distal(detections{i_im,3}(i_can))
            true_status = 'Non-distal: ';
        elseif is_undefined(detections{i_im,3}(i_can))
            true_status = 'Undefined: ';
        end
    else
        true_status = 'FP: ';
    end
    
    figure; 
    subplot(1,2,1);
    imgray(nailfold(aam_can.sr: aam_can.er, aam_can.sc:aam_can.ec));
    plot(aam_can.vessel_xy(:,1), aam_can.vessel_xy(:,2), 'g-x');
    plot(aam_can.fitted_vessel_xy(:,1), aam_can.fitted_vessel_xy(:,2), 'r-x');
    title([pred_status num2str(aam_can.model_score)]);
    
    subplot(1,2,2);
    imgray(vessel_prob(aam_can.sr: aam_can.er, aam_can.sc:aam_can.ec));
    plot(aam_can.vessel_xy(:,1), aam_can.vessel_xy(:,2), 'g-x');
    plot(aam_can.fitted_vessel_xy(:,1), aam_can.fitted_vessel_xy(:,2), 'r-x');
    title([true_status num2str(candidate_scores(i_can))]);
end
%%
width_dir = 'C:\isbe\nailfold\data\rsa_study\test\predictions\width\rf_regression\257847\';
apex_gt_dir = 'C:\isbe\nailfold\data\rsa_study\test\apex_gt\';
clear varargin;

im_list = dir('C:\isbe\nailfold\data\rsa_study\test\images\*.mat');

all_pred_widths = zeros(0,1);
all_gt_widths = zeros(0,1);

for i_im = 1:length(im_list);
    
    im_name = im_list(i_im).name(1:6);
    
    load([apex_gt_dir im_name '_gt']);
    
    vessel_width = u_load([width_dir im_name '_pred.mat']);
    vessel_width = imfilter(vessel_width, fspecial('gaussian', [10 2]));
    
    all_pred_widths = [all_pred_widths; ...
        interp2(vessel_width, apex_xy(is_distal,1), apex_xy(is_distal, 2))]; %#ok
    
    all_gt_widths = [all_gt_widths; apex_widths(is_distal)]; %#ok
end

figure; 
plot(all_gt_widths, all_pred_widths, 'rx');
%%
pred_dir = 'C:\isbe\nailfold\data\rsa_study\set12g\predictions\width\rf_regression\257847\';
width_dir = 'C:\isbe\nailfold\data\rsa_study\set12g\width_maps\';
v_mask_dir = 'C:\isbe\nailfold\data\rsa_study\set12g\vessel_centre_masks\';
im_list = dir([pred_dir '*.mat']);

all_pred_widths = zeros(0,1);
all_gt_widths = zeros(0,1);

discard = length('_pred.mat');

for i_im = 1:length(im_list)
    im_name = im_list(i_im).name(1:end-discard);
    
    pred_width = u_load([pred_dir im_name '_pred.mat']);
    gt_width = u_load([width_dir im_name '_width.mat']);
    v_mask = u_load([v_mask_dir im_name '_v_cmask']);
    
    all_pred_widths = [all_pred_widths; pred_width(v_mask)]; %#ok    
    all_gt_widths = [all_gt_widths; gt_width(v_mask)]; %#ok
end

figure; axis equal; hold on;
plot(all_gt_widths, all_pred_widths, 'rx');

figure; 
hist(all_pred_widths - all_gt_widths, -100:100);

figure; 
plot((all_pred_widths + all_gt_widths)/2, all_pred_widths - all_gt_widths, 'rx');

counts = hist3([(all_pred_widths - all_gt_widths) (all_pred_widths + all_gt_widths)/2], {-100:100, 0:200});
figure; imgray(log(counts)); colormap(jet);
%%
missed = [];
im_list = dir('C:\isbe\nailfold\data\rsa_study\test\images\*.mat');
candidates_dir = 'C:\isbe\nailfold\data\rsa_study\test\apex_maps\frog\post_merged\';
for i_im = 1:length(im_list);
    kept = NaN;
    im_name = im_list(i_im).name(1:6);
    load([candidates_dir im_name '_candidates'], 'kept');
    if isnan(kept)
        missed(end+1) = i_im;
    end
end
