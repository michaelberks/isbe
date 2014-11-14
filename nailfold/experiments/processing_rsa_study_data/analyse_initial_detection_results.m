%_------------------------------------------------------------------------
% Script for analysing the reults of detecting apices in the test set of
% images
%-------------------------------------------------------------------------
%%
%-------------------------------------------------------------------------
%% 1) First make sure we have ground truth position for apices marked

rsa_dir = 'rsa_study/';

centre_dir = [nailfoldroot 'data/' rsa_dir 'test/vessel_centres/'];
apex_gt_dir = [nailfoldroot 'data/' rsa_dir 'test/apex_gt/'];
create_folder(apex_gt_dir);

markup_names = dir([nailfoldroot 'data/' rsa_dir 'markup/']);
markup_dirs = cell(0,1);
for i_dir = 3:length(markup_names)
    if markup_names(i_dir).isdir
        markup_dirs{end+1,1} = [nailfoldroot 'data/' rsa_dir 'markup/' markup_names(i_dir).name '/']; %#ok
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
        distal_size = cell(0,1);
        distal_shape = cell(0,1);
        
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
                        distal_size{end+1,:} = vessel_markup.vessels(i_v).ncm_vessel_properties.size; %#ok
                        distal_shape{end+1,:} = vessel_markup.vessels(i_v).ncm_vessel_properties.shape; %#ok

                    end
                end
            else
                %Save the position of the anchor in the non-distal holder
                non_distal_xy(end+1,:) = anchor_xy; %#ok
                
            end
        end
    end
    %
    vessel_markup = rmfield(vessel_markup, {'vessels', 'haemorrhages'});
    save([apex_gt_dir im_num '_gt.mat'], 'distal_xy', 'distal_width', 'non_distal_xy', 'undefined_xy',...
        'vessel_markup', 'distal_size', 'distal_shape');
end
%%
% For better ground truth, make up apex clusters using all available
% markers

rsa_dir = 'rsa_study/';

image_dir = [nailfoldroot 'data/' rsa_dir 'test/images/'];
markup_dir = [nailfoldroot 'data/' rsa_dir 'markup/'];
cluster_dir = [nailfoldroot 'data/' rsa_dir 'apex_clusters/'];
create_folder(cluster_dir)

markers = dir(markup_dir);
markers = markers(3:end);
discard_markers = struct2cell(markers);
markers(~cell2mat(discard_markers(4,:))) = [];
    
test_list = dir([image_dir '*.mat']);
%
for i_im = 1:length(test_list)
    
    im_num = test_list(i_im).name(1:6);      
    [vessels] = cluster_vessel_apices(im_num, markup_dir, markers, 20, 1);
    save([cluster_dir im_num '_apex_clusters.mat'], 'vessels');
end
%%
markers_per_image = zeros(length(test_list),1);
for i_im = 1:length(test_list)
    
    im_num = test_list(i_im).name(1:6);  
    load([cluster_dir im_num '_apex_clusters.mat'], 'vessels');
    markers_per_image(i_im) = length(vessels.markers);
end
%-------------------------------------------------------------------------
%% 2) Now we can run through all the images, work which apices were detected
% We also analyse:
%   - What measure detected each apex
%   - What % of detections were true
%   - For apices we missed, why were they discarded
rsa_dir = 'rsa_study/';
corr_scores_dir = [nailfoldroot 'data/' rsa_dir 'test/corr_scores/'];
shape_scores_dir = [nailfoldroot 'data/' rsa_dir 'test/shape_scores/'];
apex_gt_dir = [nailfoldroot 'data/' rsa_dir 'test/apex_gt/'];
centre_dir = [nailfoldroot 'data/' rsa_dir 'test/vessel_centres/'];
results_dir = [nailfoldroot 'data/' rsa_dir 'test/results/'];

test_list = dir([corr_scores_dir '*_cs.mat']);
num_images = length(test_list);
%%
%Set flags
do_corr = 1;
do_shape = 1;
do_all = 1;

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
    
    or_detected_distal = 0;
    or_detected_non_distal = 0;
    or_detected_undefined = 0;
    
    and_detected_distal = 0;
    and_detected_non_distal = 0;
    and_detected_undefined = 0;
end

if do_all
    all_detections = cell(num_images, 4);
    
    all_detected_distal = 0;
    all_detected_non_distal = 0;
    all_detected_undefined = 0;
  
    total_all_candidates = 0;
    correct_all_candidates = 0;
end

non_detections = cell(num_images, 3);

%
for i_im = 1:num_images;
    
    im_num = test_list(i_im).name(1:6);
    display(['Processing image ' num2str(i_im) ' of ' num2str(num_images)]);
    
    
    load([apex_gt_dir im_num '_gt.mat']);
    load([centre_dir im_num '_vc.mat']);
    total_distal = total_distal + size(distal_xy,1);
    total_non_distal = total_non_distal + size(non_distal_xy,1);
    total_undefined = total_undefined + size(undefined_xy,1);
    
    %---------------------------------------------
    if do_corr
        load([corr_scores_dir im_num '_cs.mat']);
        [candidate_label_distal corr_detections{i_im,1}] = evaluate_apex_candidates(distal_xy, maxima_corr_xy, distal_width);
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

        %maxima_shape_xy = maxima_shape_xy(maxima_shape_vals > 100,:);
        maxima_shape_xy(101:end,:) = [];
        
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
        or_detected_distal = or_detected_distal + sum(shape_detections{i_im,1} | corr_detections{i_im,1});
        or_detected_non_distal = or_detected_non_distal + sum(shape_detections{i_im,2} | corr_detections{i_im,2});
        or_detected_undefined = or_detected_undefined + sum(shape_detections{i_im,3} | corr_detections{i_im,3});
        
        and_detected_distal = and_detected_distal + sum(shape_detections{i_im,1} & corr_detections{i_im,1});
        and_detected_non_distal = and_detected_non_distal + sum(shape_detections{i_im,2} & corr_detections{i_im,2});
        and_detected_undefined = and_detected_undefined + sum(shape_detections{i_im,3} & corr_detections{i_im,3});
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
%
save([results_dir 'results_' datestr(now,30) '.mat'],...
    '*_detected_*', 'total_*', 'correct_*_candidates', '*_detections');
%%
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
display(['% distal apices detected by correlation or shape: ' num2str(or_detected_distal / total_distal)]);
display(['% non-distal apices detected by correlation or shape: ' num2str(or_detected_non_distal / total_non_distal)]);
display(['% undefined apices detected by correlation or shape: ' num2str(or_detected_undefined / total_undefined)]);
display('***');
display(['% distal apices detected by correlation and shape: ' num2str(and_detected_distal / total_distal)]);
display(['% non-distal apices detected by correlation and shape: ' num2str(and_detected_non_distal / total_non_distal)]);
display(['% undefined apices detected by correlation and shape: ' num2str(and_detected_undefined / total_undefined)]);
display('***');
display(['% distal apices detected by centre-lines: ' num2str(all_detected_distal / total_distal)]);
display(['% non-distal apices detected by centre-lines: ' num2str(all_detected_non_distal / total_non_distal)]);
display(['% undefined apices detected by centre-lines: ' num2str(all_detected_undefined / total_undefined)]);
display('*-------------------------------------------------------------*');

%--------------------------------------------------------------------------
%% 3) Analyse why non detected points were discarded
%Make histogram of non detection reasons
total_undetected = total_distal - or_detected_distal;
discard_reasons = zeros(total_undetected, 1);
undetected_counts = zeros(num_images,1);
all_discard_counts = zeros(1, 8);
all_counts = zeros(num_images,1);
curr_idx = 0;
for i_im = 1:num_images
    im_num = test_list(i_im).name(1:6);
    
    if ~isempty(non_detections{i_im,1})
        undetected_counts(i_im) = length(non_detections{i_im,1});
        idx = curr_idx + (1:undetected_counts(i_im));
        curr_idx = idx(end);
        discard_reasons(idx) = non_detections{i_im,1};
    end
    all_counts(i_im) = length(shape_detections{i_im,1});
    
    load([centre_dir im_num '_vc.mat'], 'vessel_centre_discards');   
    all_discard_counts = all_discard_counts + hist(vessel_centre_discards, 1:8);
end
discard_counts = hist(discard_reasons, 1:8);

figure; bar(1:8, discard_counts);
title('Total apexes discarded');
%set(gca, 'xticklabel', {'Edge of mosaic', 'Straight horizontal', 'Local to XC maxima',...
%    'Outside estimated distal range', 'Too low prob', 'Too low curv', 'Local to XC maxima', 'No centreline'});
figure; bar(1:8, all_discard_counts);
title('Total points discarded');

figure; bar(1:8, discard_counts ./ all_discard_counts);
title('Apexes discarded as a percentage of points discarded');

%--------------------------------------------------------------------------
%% 4) Analyse the widths of detected/non detected vessels
total_undetected = total_distal - or_detected_distal;
detected_widths = zeros(or_detected_distal,1);
missed_widths = zeros(total_undetected, 1);

curr_idx_d = 0;
curr_idx_m = 0;
for i_im = 1:num_images
    im_num = test_list(i_im).name(1:6);
    load([apex_gt_dir im_num '_gt.mat'], 'distal_width');
    
    if (image_grade_idx(i_im) == 2) || (image_grade_idx(i_im) == 7)
        continue;
    end
    
    detected = shape_detections{i_im,1} | corr_detections{i_im,1};
    
    if any(detected)
        idx = curr_idx_d + (1:sum(detected));
        curr_idx_d = idx(end);
        detected_widths(idx) = distal_width(detected);
    end
    if any(~detected)
        idx = curr_idx_m + (1:sum(~detected));
        curr_idx_m = idx(end);
        missed_widths(idx) = distal_width(~detected);
    end
end
%
width_bins = 1:2:200;
detected_bin_counts = hist(detected_widths, width_bins);
missed_bin_counts = hist(missed_widths, width_bins);

figure; bar(width_bins, [detected_bin_counts/sum(detected_bin_counts); missed_bin_counts/sum(missed_bin_counts)]'); hold on;
plot(width_bins, missed_bin_counts/sum(missed_bin_counts), 'r', 'linewidth', 2);
plot(width_bins, detected_bin_counts/sum(detected_bin_counts), 'b', 'linewidth', 2);
title('Distribution of widths for detected and missed apices');

figure; plot(width_bins,  detected_bin_counts./ (detected_bin_counts+missed_bin_counts));
title('Detection rate as a function of vessel width');
ylabel('% of vessels detected');
xlabel('Vessel width');

detection_ratios = zeros(200,1);

for ii = 2:200
    detected_i = sum(detected_widths < ii);
    missed_i = sum(missed_widths < ii);
    detection_ratios(ii) = detected_i / (detected_i + missed_i);    
end
figure; hold on;
plot(1:200, detection_ratios, '-');
title('Detection rate of all vessels up to a given width');
ylabel('% of vessels detected');
xlabel('Vessel width');
%--------------------------------------------------------------------------
%% Analyse detections as a function of image type
image_grade = cell(num_images,1); 
marker = cell(num_images,1);
for i_im = 1:num_images
    im_num = test_list(i_im).name(1:6);
    
    load([apex_gt_dir im_num '_gt.mat'], 'vessel_markup');
    image_grade{i_im} = vessel_markup.image_grade;
    marker{i_im} = vessel_markup.observer;
end

[image_grade_idx image_grade_labels] = grp2idx(image_grade);
image_grade_counts = hist(image_grade_idx, 1:7);
[marker_idx marker_labels] = grp2idx(marker);
marker_counts = hist(marker_idx, 1:3);

%
detection_by_grade = zeros(1,7);
vessels_by_grade = zeros(1,7);

for i_im = 1:num_images
   detection_by_grade(image_grade_idx(i_im)) =  ...
       detection_by_grade(image_grade_idx(i_im)) + sum(shape_detections{i_im,1} | corr_detections{i_im,1});
   
   vessels_by_grade(image_grade_idx(i_im)) =  ...
       vessels_by_grade(image_grade_idx(i_im)) + length(shape_detections{i_im});
end
   
for i_gr = 1:7
    display(image_grade_labels{i_gr});
    display(['% distal apices detected by correlation or shape: ' num2str(detection_by_grade(i_gr) / vessels_by_grade(i_gr))]);
end
%--------------------------------------------------------------------------
%% Analyse detections by vessel size and shape
shape_labels = {'Non-specific', 'Normal', 'Angiogenic','Meandering'};
size_labels = {'Normal', 'Enlarged', 'Giant', 'Irregular', 'Undefined'};
detections_by_size = zeros(1,5);
detections_by_shape = zeros(1,4);
counts_by_size = zeros(1,5);
counts_by_shape = zeros(1,4);

for i_im = 1:num_images
    im_num = test_list(i_im).name(1:6);
    
    load([apex_gt_dir im_num '_gt.mat'], 'distal_size', 'distal_shape', 'vessel_markup');
    detected = shape_detections{i_im,1} | corr_detections{i_im,1};
    
    if (image_grade_idx(i_im) == 2) || (image_grade_idx(i_im) == 7)
        continue;
    end
    
    for i_ap = 1:length(distal_size)
        i_sz = find(strcmp(distal_size{i_ap}, size_labels));
        detections_by_size(i_sz) = detections_by_size(i_sz) + detected(i_ap);
        counts_by_size(i_sz) = counts_by_size(i_sz) + 1;
        
        i_sh = find(strcmp(distal_shape{i_ap}, shape_labels));
        detections_by_shape(i_sh) = detections_by_shape(i_sh) + detected(i_ap);
        counts_by_shape(i_sh) = counts_by_shape(i_sh) + 1;
    end
end
        
for i_sz = 1:5
    display(size_labels{i_sz});
    display(['% distal apices detected by correlation or shape: ' num2str(detections_by_size(i_sz) / counts_by_size(i_sz))]);
end   

for i_sh = 1:4
    display(shape_labels{i_sh});
    display(['% distal apices detected by correlation or shape: ' num2str(detections_by_shape(i_sh) / counts_by_shape(i_sh))]);
end
    

%--------------------------------------------------------------------------
%% 4) Save a region around all the undetected apexes
image_dir = [nailfoldroot 'data/' rsa_dir 'test/images/'];
prob_dir = [nailfoldroot 'data/' rsa_dir 'test/predictions/detection/rf_classification/222836/'];
ori_dir = [nailfoldroot 'data/' rsa_dir 'test/predictions/orientation/rf_regression/222835/'];
missed_apex_dir = [nailfoldroot 'data/' rsa_dir 'test/results/missed_apexes/'];
create_folder([missed_apex_dir 'images/']);
create_folder([missed_apex_dir 'vessel_prob/']);
create_folder([missed_apex_dir 'vessel_ori/']);
create_folder([missed_apex_dir 'vessel_centres/']);

for i_im = 1:num_images
    im_num = test_list(i_im).name(1:6);
    
    num_missed = length(non_detections{i_im,1});
    if num_missed
        load([apex_gt_dir im_num '_gt.mat'], 'distal_xy');
        load([centre_dir im_num '_vc.mat']); 
        load([shape_scores_dir im_num '_ss.mat']);
        load([corr_scores_dir im_num '_cs.mat']);
        
        apex_prior_scores_sum = squeeze(sum(sum(apex_prior_scores,2)));
        
        missed_xy = distal_xy(~(shape_detections{i_im,1} | corr_detections{i_im,1}),:);
        
        for i_ap = 1:num_missed
            
            ma_x = round(missed_xy(i_ap,1));
            ma_y = round(missed_xy(i_ap,2));
            
            %Get centre points in the region
            include_pts = ...
                (abs(vessel_centre_x - ma_x) < 100) & ...
                (abs(vessel_centre_y - ma_y) < 100);
            
            patch_centre_x = vessel_centre_x(include_pts) - ma_x + 101;
            patch_centre_y = vessel_centre_y(include_pts) - ma_y + 101;
            
            %Get shape points from the region
            include_pts = ...
                (abs(apex_x - ma_x) < 100) & ...
                (abs(apex_y - ma_y) < 100);
            
            patch_centre_shape_x = apex_x(include_pts) - ma_x + 101;
            patch_centre_shape_y = apex_y(include_pts) - ma_y + 101;
            patch_centre_shape_sum = apex_prior_scores_sum(include_pts);
            
            %Get shape maxima in the region
            include_shape_maxima = ...
                (abs(maxima_shape_xy(:,1) - ma_x) < 100) & ...
                (abs(maxima_shape_xy(:,2) - ma_y) < 100);
            
            patch_shape_maxima_x = maxima_shape_xy(include_shape_maxima,1) - ma_x + 101;
            patch_shape_maxima_y = maxima_shape_xy(include_shape_maxima,2) - ma_y + 101;
            
            %Get correlation maxima in the region
            include_corr_maxima = ...
                (abs(maxima_corr_xy(:,1) - ma_x) < 100) & ...
                (abs(maxima_corr_xy(:,2) - ma_y) < 100);
            
            patch_corr_maxima_x = maxima_shape_xy(include_corr_maxima,1) - ma_x + 101;
            patch_corr_maxima_y = maxima_shape_xy(include_corr_maxima,2) - ma_y + 101;
            
            %Take patches from the nailfold image, vessel probability and
            %orientation prediction maps
            patch_centre_prob = vessel_centre_prob(include_pts);
            patch_centre_ori = vessel_centre_ori(include_pts);
            patch_centre_curv = vessel_centre_curv(include_pts);
            patch_centre_corr = vessel_centre_corr(include_pts);
            patch_centre_discards = vessel_centre_discards(include_pts);
            
            patch_im = sample_window(vessel_im, 201, ma_y, ma_x);
            patch_prob = sample_window(vessel_prob, 201, ma_y, ma_x);
            patch_ori = sample_window(vessel_ori, 201, ma_y, ma_x);
            
            imwrite(uint8(patch_im), [missed_apex_dir 'images/' im_num '_' zerostr(i_ap,2) '.png']);
            save([missed_apex_dir 'vessel_prob/' im_num '_' zerostr(i_ap,2) '_prob.mat'], 'patch_prob');
            save([missed_apex_dir 'vessel_ori/' im_num '_' zerostr(i_ap,2) '_ori.mat'], 'patch_ori');
            save([missed_apex_dir 'vessel_centres/' im_num '_' zerostr(i_ap,2) '_vc.mat'], 'patch_centre_*', 'patch_*_maxima_*');
        end
    end
end
%%
missed_list = dir([missed_apex_dir 'images/*.png']);
dot_color = hsv(8);
r_idx = randperm(1000);
r_idx = r_idx(1:20) + 20;
for i_im = r_idx
    
    patch_im = imread([missed_apex_dir 'images/' missed_list(i_im).name]);
    patch_prob = u_load([missed_apex_dir 'vessel_prob/' missed_list(i_im).name(1:end-4) '_prob.mat']);
    load([missed_apex_dir 'vessel_centres/' missed_list(i_im).name(1:end-4) '_vc.mat']);
    
    figure;
    ax(1) = subplot(1,2,1);
    imgray(patch_im);
    
    ax(2) = subplot(1,2,2);
    imgray(patch_prob);
    
    for i_di = 0:7
        include_pts = patch_centre_discards == i_di;
        %plot(ax(1), patch_centre_x(include_pts), patch_centre_y(include_pts), '.', 'markersize', 6, 'markeredgecolor', dot_color(i_di+1,:));
        plot(ax(2), patch_centre_x(include_pts), patch_centre_y(include_pts), '.', 'markersize', 6, 'markeredgecolor', dot_color(i_di+1,:));
    end
    
    for i_ax = 1:2
        plot(ax(i_ax), 101, 101, 'c*', 'markersize', 12);
    end
    plot(ax(2), patch_corr_maxima_x, patch_corr_maxima_y, 'gs', 'markersize', 12); 
    plot(ax(2),patch_shape_maxima_x, patch_shape_maxima_y, 'go', 'markersize', 12);
        
end
%%
shape_prior_bins = 1:5:400;
shape_prior_counts = zeros(2,length(shape_prior_bins));
for i_im = 1:num_images
    
    im_num = test_list(i_im).name(1:6);
    
    load([shape_scores_dir im_num '_ss.mat']);
    shape_prior_counts(1,:) = shape_prior_counts(1,:) + ...
        hist(maxima_shape_vals(shape_detections{i_im,4}), shape_prior_bins);
    
    shape_prior_counts(2,:) = shape_prior_counts(2,:) + ...
        hist(maxima_shape_vals(~shape_detections{i_im,4}), shape_prior_bins);
    
end

figure; 
subplot(1,2,1); bar(shape_prior_bins, shape_prior_counts(1,:));
subplot(1,2,2); bar(shape_prior_bins, shape_prior_counts(2,:));

            
            
        
        
        
        
        




