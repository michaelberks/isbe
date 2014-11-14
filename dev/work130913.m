rsa_dir = 'rsa_study/';
corr_pos_dir = [nailfoldroot 'data/' rsa_dir 'test/corr_pos/'];
test_list = dir([corr_pos_dir '*_cp.mat']);
%
num_images = length(test_list);

for i_im = 1:20
    
    im_num = test_list(i_im).name(1:6);
    
    load([corr_pos_dir im_num '_cp.mat']);

    %all_xy = [distal_xy; undefined_xy; non_distal_xy];
    
    all_xy = maxima_corr_pos;

    num_candidates = size(all_xy,1);
    
    if num_candidates < 10
        continue;
    end

    violations_mat = true(num_candidates, num_candidates);
    tan_mat = zeros(num_candidates, num_candidates);
    for i_can = 1:num_candidates

        leave_out_idx = setdiff(1:num_candidates, i_can);

        tan_vectors = abs(atan(...
            (all_xy(leave_out_idx,2)-all_xy(i_can,2)) ./...
            (all_xy(leave_out_idx,1)-all_xy(i_can,1)) ));

        violations_mat(leave_out_idx, i_can) = tan_vectors < pi/4;
        tan_mat(leave_out_idx,i_can) = tan_vectors;
        
    end

    happy_list = 1:num_candidates;
    violations_i = violations_mat;
    tan_i = tan_mat;
    while any(~violations_i(:))
        %vio_count = sum(violations_i);
        %[~, worst_offender] = min(vio_count);
        tan_sum = sum(tan_mat);
        [~, worst_offender] = max(tan_mat);
        
        violations_i(worst_offender,:) = [];
        violations_i(:,worst_offender) = [];
        
        tan_mat(worst_offender,:) = [];
        tan_mat(:,worst_offender) = [];
        
        happy_list(worst_offender) = [];
        
    end

    %figure;
    %subplot(1,2,1); imgray(violations_mat);
    %subplot(1,2,2); imgray(violations_i);

    figure; hold all; axis equal ij;
    plot(distal_xy(:,1), distal_xy(:,2), 'r^'); 
    plot(non_distal_xy(:,1), non_distal_xy(:,2), 'bs'); 
    plot(undefined_xy(:,1), undefined_xy(:,2), 'r*'); 
    plot(all_xy(happy_list,1), all_xy(happy_list,2), 'go');
    plot(maxima_corr_pos(:,1), maxima_corr_pos(:,2), 'c+');
end
%%
rsa_dir = 'rsa_study/';

centre_dir = [nailfoldroot 'data/' rsa_dir 'test/vessel_centres/'];
apex_gt_dir = [nailfoldroot 'data/' rsa_dir 'test/apex_gt/'];
create_folder(apex_gt_dir);

markup_names = dir([nailfoldroot 'data/' rsa_dir 'markup/']);
markup_dirs = cell(0,1);
for i_dir = 3:length(markup_names)
    if markup_names(i_dir).isdir
        markup_dirs{end+1,1} = [nailfoldroot 'data/' rsa_dir 'markup/' markup_names(i_dir).name '/'];
    end
end
%
test_list = dir([centre_dir '*_vc.mat']);
%
for i_im = 1:length(test_list)
    
    im_num = test_list(i_im).name(1:6);
    
    for i_dir = 1:length(markup_dirs)
        vessel_markup_list_i = dir([markup_dirs{i_dir} '*' im_num '*.txt']);
        if ~isempty(vessel_markup_list_i)
            markup_dir = markup_dirs{i_dir};
        	vessel_markup_list = vessel_markup_list_i;
        end
    end
    
    if ~isempty(vessel_markup_list)
        
        distal_xy = zeros(0,2);
        distal_width = zeros(0,1);
        non_distal_xy = zeros(0,2);
        undefined_xy = zeros(0,2);
        
        %Read in markup
        vessel_markup = read_markup_from([markup_dir vessel_markup_list(end).name]);
        
        %Loop through each marked vessel apex
        num_vessels = length(vessel_markup.vessels);
        for i_v = 1:num_vessels

            %Check this is a valid vessel
            anchor_xy = vessel_markup.vessels(i_v).anchor;

            if isempty(anchor_xy); continue; end

            %Check if distal
            is_distal = vessel_markup.vessels(i_v).ncm_vessel_properties.is_distal;

            if is_distal
                
                num_apices = length(vessel_markup.vessels(i_v).apices);
                for i_a = 1:num_apices
                    if isempty(vessel_markup.vessels(i_v).apices(i_a).inner_point)                    
                        %Save the position of the anchor in the undefined holder
                        undefined_xy(end+1,:) = anchor_xy; %#ok

                    else
                        %Compute the centre of the apex and save in distal
                        %holder
                        apex_xy =  ...
                            [ vessel_markup.vessels(i_v).apices(i_a).outer_point;...
                              vessel_markup.vessels(i_v).apices(i_a).inner_point];
                        apex_width = sqrt(sum(diff(apex_xy).^2));
                        
                        distal_xy(end+1,:) = mean(apex_xy); %#ok 
                        distal_width(end+1,:) = apex_width; %#ok

                    end
                end
            else
                %Save the position of the anchor in the non-distal holder
                non_distal_xy(end+1,:) = anchor_xy; %#ok
                
            end
        end
    end
    %
    save([apex_gt_dir im_num '_gt.mat'], 'distal_xy', 'distal_width', 'non_distal_xy', 'undefined_xy');
end
%%
rsa_dir = 'rsa_study/';
corr_scores_dir = [nailfoldroot 'data/' rsa_dir 'test/corr_scores/'];
shape_scores_dir = [nailfoldroot 'data/' rsa_dir 'test/shape_scores/'];
apex_gt_dir = [nailfoldroot 'data/' rsa_dir 'test/apex_gt/'];
centre_dir = [nailfoldroot 'data/' rsa_dir 'test/vessel_centres/'];

test_list = dir([corr_scores_dir '*_cs.mat']);
num_images = length(test_list);
do_corr = 0;
do_shape = 0;
do_all = 0;

total_distal = 0;
total_non_distal = 0;
total_undefined = 0;

if do_corr
    corr_detections = cell(num_images, 4);
    
    corr_detected_distal = 0;
    corr_detected_non_distal = 0;
    corr_detected_undefined = 0;
    
    total_corr_candidates = 0;
    correct_corr_candidates = 0;
end

if do_shape
    shape_detections = cell(num_images, 4);
    
    shape_detected_distal = 0;
    shape_detected_non_distal = 0;
    shape_detected_undefined = 0;
    
    total_shape_candidates = 0;
    correct_shape_candidates = 0;
end

if do_corr && do_shape
    all_detections = cell(num_images, 4);
    
    detected_distal = 0;
    detected_non_distal = 0;
    detected_undefined = 0;
end

if do_all
    all_detected_distal = 0;
    all_detected_non_distal = 0;
    all_detected_undefined = 0;
  
    total_all_candidates = 0;
    correct_all_candidates = 0;
end

non_detections = cell(num_images, 4);


for i_im = 1:320%num_images;
    
    im_num = test_list(i_im).name(1:6);
    display(['Processing image ' num2str(i_im) ' of ' num2str(num_images)]);
    
    
    load([apex_gt_dir im_num '_gt.mat'], 'distal_xy', 'distal_width', 'non_distal_xy', 'undefined_xy');
    load([centre_dir im_num '_vc.mat']);
    total_distal = total_distal + size(distal_xy,1);
    total_non_distal = total_non_distal + size(non_distal_xy,1);
    total_undefined = total_undefined + size(undefined_xy,1);
    
    %---------------------------------------------
    if do_corr
        load([corr_scores_dir im_num '_cs.mat']);
        [candidate_label_distal corr_detections{i_im,1}] = evaluate_apex_candidates(distal_xy, maxima_corr_xy, 15);
        [candidate_label_non_distal corr_detections{i_im,2}] = evaluate_apex_candidates(non_distal_xy, maxima_corr_xy, 15);
        [candidate_label_undefined corr_detections{i_im,3}] = evaluate_apex_candidates(undefined_xy, maxima_corr_xy, 15);

        corr_detections{i_im,4} = ...
            candidate_label_distal | ...
            candidate_label_non_distal | ...
            candidate_label_undefined;

        corr_detected_distal = corr_detected_distal + sum(corr_detections{i_im,1});
        corr_detected_non_distal = corr_detected_non_distal + sum(corr_detections{i_im,2});
        corr_detected_undefined = corr_detected_undefined + sum(corr_detections{i_im,3});

        total_corr_candidates = total_corr_candidates + length(corr_detections{i_im,4});
        correct_corr_candidates = correct_corr_candidates + sum(corr_detections{i_im,4});
    end
    %-------------------------------------------------------
    if do_shape
        load([shape_scores_dir im_num '_ss.mat']);

        [candidate_label_distal shape_detections{i_im,1}] = evaluate_apex_candidates(distal_xy, maxima_shape_xy, 15);
        [candidate_label_non_distal shape_detections{i_im,2}] = evaluate_apex_candidates(non_distal_xy, maxima_shape_xy, 15);
        [candidate_label_undefined shape_detections{i_im,3}] = evaluate_apex_candidates(undefined_xy, maxima_shape_xy, 15);

        shape_detections{i_im,4} = ...
            candidate_label_distal | ...
            candidate_label_non_distal | ...
            candidate_label_undefined;

        shape_detected_distal = shape_detected_distal + sum(shape_detections{i_im,1});
        shape_detected_non_distal = shape_detected_non_distal + sum(shape_detections{i_im,2});
        shape_detected_undefined = shape_detected_undefined + sum(shape_detections{i_im,3});

        total_shape_candidates = total_shape_candidates + length(shape_detections{i_im,4});
        correct_shape_candidates = correct_shape_candidates + sum(shape_detections{i_im,4});
    end
    
    %----------------------------------------------------------
    
    if do_corr && do_shape
        detected_distal = detected_distal + sum(shape_detections{i_im,1} | corr_detections{i_im,1});
        detected_non_distal = detected_non_distal + sum(shape_detections{i_im,2} | corr_detections{i_im,2});
        detected_undefined = detected_undefined + sum(shape_detections{i_im,3} | corr_detections{i_im,3});
    end
    %----------------------------------------------------------
    if do_all

        [candidate_label_distal all_detections{i_im,1}] = evaluate_apex_candidates(distal_xy, [vessel_centre_x vessel_centre_y], 15);
        [candidate_label_non_distal all_detections{i_im,2}] = evaluate_apex_candidates(non_distal_xy, [vessel_centre_x vessel_centre_y], 15);
        [candidate_label_undefined all_detections{i_im,3}] = evaluate_apex_candidates(undefined_xy, [vessel_centre_x vessel_centre_y], 15);

        all_detected_distal = all_detected_distal + sum(all_detections{i_im,1});
        all_detected_non_distal = all_detected_non_distal + sum(all_detections{i_im,2});
        all_detected_undefined = all_detected_undefined + sum(all_detections{i_im,3});

        total_all_candidates = total_all_candidates + length(all_detections{i_im,4});
        correct_all_candidates = correct_all_candidates + sum(all_detections{i_im,4});
    end
    %--------------------------------------------------------------------
    undetected = ~(shape_detections{i_im,1} | corr_detections{i_im,1});
    
    undetected_x = distal_xy(undetected,1);
    undetected_y = distal_xy(undetected,2);
    num_pts = length(undetected_x);
    non_detections{i_im,1} = zeros(num_pts,1);
    for i_pt = 1:num_pts
        %find the nearest centreline pixel
        dists = (vessel_centre_x - undetected_x(i_pt)).^2 + (vessel_centre_y - undetected_y(i_pt)).^2;
        [min_dist, min_idx] = min(dists);
        if sqrt(min_dist) <= 15
            non_detections{i_im,1}(i_pt) = vessel_centre_discards(min_idx);
        else
            non_detections{i_im,1}(i_pt) = 8;
        end
    end
    non_detections{i_im,4} = distal_width(undetected);
    %
    undetected = ~(shape_detections{i_im,2} | corr_detections{i_im,2});
    
    undetected_x = non_distal_xy(undetected,1);
    undetected_y = non_distal_xy(undetected,2);
    num_pts = length(undetected_x);
    non_detections{i_im,2} = zeros(num_pts,1);
    for i_pt = 1:num_pts
        %find the nearest centreline pixel
        dists = (vessel_centre_x - undetected_x(i_pt)).^2 + (vessel_centre_y - undetected_y(i_pt)).^2;
        [min_dist, min_idx] = min(dists);
        if sqrt(min_dist) <= 15
            non_detections{i_im,2}(i_pt) = vessel_centre_discards(min_idx);
        end
    end
    %
    undetected = ~(shape_detections{i_im,3} | corr_detections{i_im,3});
    
    undetected_x = undefined_xy(undetected,1);
    undetected_y = undefined_xy(undetected,2);
    num_pts = length(undetected_x);
    non_detections{i_im,3} = zeros(num_pts,1);
    for i_pt = 1:num_pts
        %find the nearest centreline pixel
        dists = (vessel_centre_x - undetected_x(i_pt)).^2 + (vessel_centre_y - undetected_y(i_pt)).^2;
        [min_dist, min_idx] = min(dists);
        if sqrt(min_dist) <= 15
            non_detections{i_im,3}(i_pt) = vessel_centre_discards(min_idx);
        end
    end
    %----------------------------------------------------------------------    
    
end
display('*-------------------------------------------------------------*');    
display(['% distal apices detected by correlation: ' num2str(corr_detected_distal / total_distal)]);
display(['% non-distal apices detected by correlation: ' num2str(corr_detected_non_distal / total_non_distal)]);
display(['% undefined apices detected by correlation: ' num2str(corr_detected_undefined / total_undefined)]);
display(['% correlation maxima that were apices: ' num2str(correct_corr_candidates / total_corr_candidates)]);
display('***');
display(['% distal apices detected by shape prior scores: ' num2str(shape_detected_distal / total_distal)]);
display(['% non-distal apices detected by shape prior scores: ' num2str(shape_detected_non_distal / total_non_distal)]);
display(['% undefined apices detected by shape prior scores: ' num2str(shape_detected_undefined / total_undefined)]);
display(['% shape prior maxima that were apices: ' num2str(correct_shape_candidates / total_shape_candidates)]);
display('***');
display(['% distal apices detected by correlation or shape: ' num2str(detected_distal / total_distal)]);
display(['% non-distal apices detected by correlation or shape: ' num2str(detected_non_distal / total_non_distal)]);
display(['% undefined apices detected by correlation or shape: ' num2str(detected_undefined / total_undefined)]);
display('***');
display(['% distal apices detected by centre-lines: ' num2str(all_detected_distal / total_distal)]);
display(['% non-distal apices detected by centre-lines: ' num2str(all_detected_non_distal / total_non_distal)]);
display(['% undefined apices detected by centre-lines: ' num2str(all_detected_undefined / total_undefined)]);
display('*-------------------------------------------------------------*');

%%
%Make histogram of non detection reasons
total_undetected = total_distal - detected_distal;
discard_reasons = zeros(total_undetected, 1);
missed_widths = zeros(total_undetected, 1);
undetected_counts = zeros(num_images,1);
all_counts = zeros(num_images,1);
curr_idx = 0;
for i_im = 1:num_images
    if ~isempty(non_detections{i_im,1})
        undetected_counts(i_im) = length(non_detections{i_im,1});
        idx = curr_idx + (1:undetected_counts(i_im));
        curr_idx = idx(end);
        discard_reasons(idx) = non_detections{i_im,1};
        missed_widths(idx) = non_detections{i_im,4};
    end
    all_counts(i_im) = length(shape_detections{i_im,1});
end

figure; hist(discard_reasons, 1:8);
discard_counts = hist(discard_reasons, 1:8);
%%
all_discard_counts = zeros(1, 8);
for i_im = 1:num_images
    im_num = test_list(i_im).name(1:6);
    load([centre_dir im_num '_vc.mat'], 'vessel_centre_discards');
    
    all_discard_counts = all_discard_counts + hist(vessel_centre_discards, 1:8);
end
figure; bar(1:8, all_discard_counts);
%%

and_detected_distal = 0;
and_detected_non_distal = 0;
and_detected_undefined = 0;
for i_im = 1:num_images;
    and_detected_distal = and_detected_distal + sum(shape_detections{i_im,1} & corr_detections{i_im,1});
    and_detected_non_distal = and_detected_non_distal + sum(shape_detections{i_im,2} & corr_detections{i_im,2});
    and_detected_undefined = and_detected_undefined + sum(shape_detections{i_im,3} & corr_detections{i_im,3});
end
display(['% distal apices detected by correlation and shape: ' num2str(and_detected_distal / total_distal)]);
display(['% non-distal apices detected by correlation and shape: ' num2str(and_detected_non_distal / total_non_distal)]);
display(['% undefined apices detected by correlation and shape: ' num2str(and_detected_undefined / total_undefined)]);

%%
for i_im = 3%sorted_i(sorted_p == 1)'%find(undetected_counts > 10)'
    im_num = test_list(i_im).name(1:6);

    load([apex_gt_dir im_num '_gt.mat'], 'distal_xy', 'distal_width', 'non_distal_xy', 'undefined_xy');
    load(['C:\isbe\nailfold\data\rsa_study\test\vessel_centres\' im_num '_vc.mat'], 'vessel_centre_y', 'vessel_centre_x', 'vessel_centre_discards', 'nrows', 'ncols');
    load(['C:\isbe\nailfold\data\rsa_study\test\images\' im_num '.mat']);

    undetected = ~(shape_detections{i_im,1} | corr_detections{i_im,1});
    corr_detected = corr_detections{i_im,1};
    shape_detected = shape_detections{i_im,1};
    
    figure; imgray(nailfold);
    caxis([min(nailfold(vessel_centre_idx)) max(nailfold(vessel_centre_idx))]);
    plot(distal_xy(:,1), distal_xy(:,2), 'rx');
    plot(distal_xy(corr_detected,1), distal_xy(corr_detected,2), 'go');
    plot(distal_xy(shape_detected,1), distal_xy(shape_detected,2), 'cs');
    plot(distal_xy(undetected,1), distal_xy(undetected,2), 'rx');
end
%%
i_im = 3;
im_num = test_list(i_im).name(1:6);
load(['C:\isbe\nailfold\data\rsa_study\test\vessel_centres\' im_num '_vc.mat'], 'vessel_centre_y', 'vessel_centre_x', 'vessel_centre_discards', 'nrows', 'ncols');
vessel_centre_idx = sub2ind([nrows ncols], vessel_centre_y, vessel_centre_x);
    
discard_reason_map = zeros(nrows, ncols);   
discard_reason_map(vessel_centre_idx) = vessel_centre_discards;

custom_map = [0 0 0; hsv(6); 1 1 1];
figure; imgray(discard_reason_map); colormap(custom_map); colorbar;
plot(distal_xy(:,1), distal_xy(:,2), 'rx');
%%
do_plot = 1;
num_angles = 11;
bin_centres = -240:40:240;
ref_angles = linspace(-pi/4, pi/4, num_angles);

feature_vectors = zeros(num_images,num_angles,length(bin_centres)-2);
targets = nan(num_images, 2);

max_angles = zeros(num_images,1);
offset_means = zeros(num_images,1);
%%
for i_im = 11:20%num_images;%[ 6 40 92 202 219 223 255 266 323 344 494 517 546]%
    im_num = test_list(i_im).name(1:6);
    load([apex_gt_dir im_num '_gt.mat'], 'distal_xy', 'distal_width', 'non_distal_xy', 'undefined_xy');
    load(['C:\isbe\nailfold\data\rsa_study\test\vessel_centres\' im_num '_vc.mat']);
    
    distal_row_x = [distal_xy(:,1); undefined_xy(:,1)] - ncols/2;
    if length(distal_row_x) < 5
        continue;
    end
    
    distal_row_y = [distal_xy(:,2); undefined_xy(:,2)] - nrows/2;
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
figure;
plot(atan(targets(:,1)), max_angles, 'x');

figure;
plot(targets(:,2), offset_means, 'x');
%%
keep_pts = ~isnan(targets(:,1)) & abs(targets(:,1))<pi/4;
y = targets(keep_pts,:);
X = feature_vectors(keep_pts,:);

rf_args.prediction_type = 'rf_regression';
rf_args.n_trees = 100;
rf_args.d = [];
rf_args.w_prior = 0;
rf_args.impure_thresh = 1.0000e-008;
rf_args.split_min = 1;
rf_args.end_cut_min = 0;
rf_args.do_ubound = 0;
rf_args.quiet = 1;
rf_args.do_circular = [];
rf_args.overwrite = 1;
rf_args.minimise_size = 0;
rf_args.split_criterion = 'ssq';
rf_args.var_criterion = 'ssq';

rf_args.sampling_args.sampling_method = 'sample_saved_training_data';
rf_args.decomposition_args = [];

%Train
rf_dir = ['C:\isbe\nailfold\data\rsa_study\models\distal_row_rfs\' datestr(now, 30) '\'];
rf_args.tree_dir = [rf_dir 'trees/'];
rf_args.sampling_args.X = X(1:300,:);


rf_args.sampling_args.y = atan(y(1:300, 1));
angle_predictor = random_forest_reg_train(rf_args);
angle_p = random_forest_reg_predict(angle_predictor, X(301:end, :));

figure; plot(atan(y(301:end, 1)), angle_p, 'x');
axis equal; hold on;
plot([-1 1], [-1 1], 'r');

rf_args.sampling_args.y = y(1:300, 2);
offset_predictor = random_forest_reg_train(rf_args);
offset_p = random_forest_reg_predict(offset_predictor, X(301:end, :));
figure; plot(y(301:end, 2), offset_p, 'x');
axis equal;
%%
plot_rows = 4;
plot_cols = 6;

for i_ap = 1:67
    
    plot_num = rem(i_ap - 1, plot_cols*plot_rows) + 1;
    
    load(['C:\isbe\nailfold\data\rsa_study\test\aam\10147c\apex' zerostr(i_ap, 4) '_aam.mat']);
    %load(['C:\isbe\nailfold\data\p10598\candidates\template_matching\aligned\apex' zerostr(i_ap, 4) '_aam.mat']);
    %load(['C:\isbe\nailfold\data\d4\candidates\template_matching\aligned\apex' zerostr(i_ap, 4) '_aam.mat']);
    
    if i_ap == 1
        if strcmp(apex_candidate.image_path(end-2:end), 'mat')
            im = u_load(apex_candidate.image_path);
        else
            im = imread(apex_candidate.image_path);
            im = im(:,:,1);
        end
    end
    
    if plot_num == 1;
        figure;
    end
    
    roi = im(apex_candidate.sr:apex_candidate.er, apex_candidate.sc:apex_candidate.ec);
    subplot(plot_rows,plot_cols,plot_num); imgray(roi);
    plot(apex_candidate.vessel_xy(:,1), apex_candidate.vessel_xy(:,2), '-x');
    plot(apex_candidate.fitted_vessel_xy(:,1), apex_candidate.fitted_vessel_xy(:,2), '-x');
    title(num2str(apex_candidate.model_score));
end