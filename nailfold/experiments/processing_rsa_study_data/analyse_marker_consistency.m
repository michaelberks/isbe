%_------------------------------------------------------------------------
% Script for analysing the reults of detecting apices in the test set of
% images
%-------------------------------------------------------------------------
%%
%-------------------------------------------------------------------------
%% 1) For better ground truth, make up apex clusters using all available
% markers

rsa_dir = 'rsa_study/';

image_dir = [nailfoldroot 'data/' rsa_dir 'test/images/'];
markup_dir = [nailfoldroot 'data/' rsa_dir 'markup/'];
cluster_dir = [nailfoldroot 'data/' rsa_dir 'test/apex_clusters/'];
post_cluster_dir = [nailfoldroot 'data/' rsa_dir 'test/apex_clusters_merged/'];
apex_gt_dir = [nailfoldroot 'data/' rsa_dir 'test/apex_gt/'];
create_folder(cluster_dir);
create_folder(post_cluster_dir);

markers = dir(markup_dir);
markers = markers(3:end);
discard_markers = struct2cell(markers);
markers(~cell2mat(discard_markers(4,:))) = [];

grades = {'Non-specific', 'Normal', 'Early', 'Active', 'Late', 'Ungradeable_Quality', 'Ungradeable_Condition'};
vessel_shapes = {'Normal', 'Non-specific', 'Angiogenic', 'Meandering', 'Undefined', 'Distal'};
vessel_sizes = {'Normal', 'Enlarged', 'Irregular', 'Giant', 'Undefined', 'Distal'};
    
test_list = dir([image_dir '*.mat']);
%%
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
%%
prob_dir = [nailfoldroot 'data/' rsa_dir 'test/predictions/detection/rf_classification/257273/'];
dist_thresh = 100;
patch_sz2 = 50;
%
total_to_correct = 0;
plot_num = 1;
for i_im = 1:length(test_list);
    
    if markers_per_image(i_im) > 1
        im_num = test_list(i_im).name(1:6);  
        vessel_prob = u_load([prob_dir im_num '_pred.mat']);
        load([cluster_dir im_num '_apex_clusters.mat'], 'vessels');
        
        [vessels] = post_merge_apexes(vessels, vessel_prob, dist_thresh, patch_sz2,...
            0.5, 20, 0);
        
        %save([post_cluster_dir im_num '_apex_clusters.mat'], 'vessels');
    end
end
%%
for i_im = 1:length(test_list);
    
    im_num = test_list(i_im).name(1:6);  
    
    if exist([post_cluster_dir im_num '_apex_clusters.mat'], 'file')
        load([post_cluster_dir im_num '_apex_clusters.mat'], 'vessels');
    else
        load([cluster_dir im_num '_apex_clusters.mat'], 'vessels');
    end
    
    grades = vessels.grades;
    majority_grade = vessels.majority_grade;
    gradeable = ~any(...
        strcmpi(grades, 'Ungradeable_Quality') |...
        strcmpi(grades, 'Ungradeable_Condition') );
    
    apex_xy = vessels.cluster_centres;
    apex_widths = vessels.cluster_radius;
    num_im_markers = length(vessels.markers);
    num_apex_markers = vessels.num_markers;
    is_non_distal = strcmpi(vessels.majority_shapes, 'Distal');
    is_undefined = strcmpi(vessels.majority_shapes, 'Undefined');
    is_distal = ~is_non_distal & ~is_undefined;

    apex_shape = vessels.majority_shapes;
    apex_size = vessels.majority_sizes;
       
    %Save a new apex GT structure
    save([apex_gt_dir im_num '_gt.mat'], 'gradeable', 'apex_xy', 'apex_widths', 'num_im_markers', 'num_apex_markers',...
        'is_non_distal', 'is_undefined', 'is_distal', 'grades', 'apex_size', 'apex_shape', 'majority_grade');
end
    
%-------------------------------------------------------------------------
%% 2) Analyse image grade co-occurences
grade_coocs2 = zeros(length(grades));

for i_im = find(markers_per_image == 2)'
    im_num = test_list(i_im).name(1:6);  
    load([cluster_dir im_num '_apex_clusters.mat'], 'vessels');
    
    r_idx = find(strcmpi(vessels.grades{1}, grades));
    c_idx = find(strcmpi(vessels.grades{2}, grades));
    
    if isempty(r_idx)
        display(vessels.grades{1});
    end
    if isempty(c_idx)
        display(vessels.grades{2});
    end

    grade_coocs2(r_idx, c_idx) = grade_coocs2(r_idx, c_idx) + 1;
end
%
for i_g1 = 1:7
    for i_g2 = i_g1+1:7
        grade_coocs2(i_g1, i_g2) = grade_coocs2(i_g1, i_g2) + grade_coocs2(i_g2, i_g1);
        grade_coocs2(i_g2, i_g1) = grade_coocs2(i_g1, i_g2);
    end
end
display(grade_coocs2);
%%
grade_coocs3 = zeros(length(grades));

for i_im = find(markers_per_image == 3)'
    im_num = test_list(i_im).name(1:6);  
    load([cluster_dir im_num '_apex_clusters.mat'], 'vessels');
    
    c_idx = find(strcmpi(vessels.majority_grade, grades));
    for i_g = 1:3
        r_idx = strcmpi(vessels.grades{i_g}, grades);
        grade_coocs3(r_idx, c_idx) = grade_coocs3(r_idx, c_idx) + 1;
    end
end
display(grade_coocs3);
%%
shape_coocs2 = zeros(length(vessel_shapes));
size_coocs2 = zeros(length(vessel_sizes));

for i_im = find(markers_per_image >= 2)'
    im_num = test_list(i_im).name(1:6);  
    load([cluster_dir im_num '_apex_clusters.mat'], 'vessels');
    
    if (ismember(grades{6}, vessels.grades) || ismember(grades{7}, vessels.grades))
        continue;
    end
    
    for i_ve = 1:size(vessels.cluster_centres,1)
        
        if length(vessels.cluster_shapes{i_ve}) == 2
            
            r_idx = find(strcmpi(vessels.cluster_shapes{i_ve}{1}, vessel_shapes));
            c_idx = find(strcmpi(vessels.cluster_shapes{i_ve}{2}, vessel_shapes));
            shape_coocs2(r_idx, c_idx) = shape_coocs2(r_idx, c_idx) + 1;
            
            r_idx = find(strcmpi(vessels.cluster_sizes{i_ve}{1}, vessel_sizes));
            c_idx = find(strcmpi(vessels.cluster_sizes{i_ve}{2}, vessel_sizes));
            size_coocs2(r_idx, c_idx) = size_coocs2(r_idx, c_idx) + 1;
        end
    end
    
            
end
%
for i_s1 = 1:length(vessel_shapes)
    for i_s2 = i_s1+1:length(vessel_shapes)
        shape_coocs2(i_s1, i_s2) = shape_coocs2(i_s1, i_s2) + shape_coocs2(i_s2, i_s1);
        shape_coocs2(i_s2, i_s1) = shape_coocs2(i_s1, i_s2);
    end
end
display(shape_coocs2);

for i_s1 = 1:length(vessel_sizes)
    for i_s2 = i_s1+1:length(vessel_sizes)
        size_coocs2(i_s1, i_s2) = size_coocs2(i_s1, i_s2) + shape_coocs2(i_s2, i_s1);
        size_coocs2(i_s2, i_s1) = size_coocs2(i_s1, i_s2);
    end
end
display(size_coocs2);
%%
vessel_counts2 = zeros(2,1);
for i_im = find(markers_per_image == 2)'
    im_num = test_list(i_im).name(1:6);  
    load([cluster_dir im_num '_apex_clusters.mat'], 'vessels');
    if (ismember(grades{6}, vessels.grades) || ismember(grades{7}, vessels.grades))
        continue;
    end
    
    for i_ve = 1:size(vessels.cluster_centres,1)
        
        if ~strcmpi(vessels.majority_shapes{i_ve}, 'Distal')
            num_hits = length(unique(vessels.cluster_members{i_ve}));
            vessel_counts2(num_hits) = vessel_counts2(num_hits) + 1;
        end
    end
end
display(vessel_counts2);

vessel_counts3 = zeros(3,1);
for i_im = find(markers_per_image == 3)'
    im_num = test_list(i_im).name(1:6);  
    load([cluster_dir im_num '_apex_clusters.mat'], 'vessels');
    if (ismember(grades{6}, vessels.grades) || ismember(grades{7}, vessels.grades))
        continue;
    end
    
    for i_ve = 1:size(vessels.cluster_centres,1)
        
        if ~strcmpi(vessels.majority_shapes{i_ve}, 'Distal')
            num_hits = length(unique(vessels.cluster_members{i_ve}));
            vessel_counts3(num_hits) = vessel_counts3(num_hits) + 1;
        end
    end
end
display(vessel_counts3);
%%


vessel_counts2_shape = zeros(2,4);
vessel_counts2_size = zeros(2,5);
vessel_counts2_grade = zeros(2,7);

for i_im = find(markers_per_image == 2)'
    im_num = test_list(i_im).name(1:6);  
    load([cluster_dir im_num '_apex_clusters.mat'], 'vessels');
    if (ismember(grades{6}, vessels.grades) || ismember(grades{7}, vessels.grades))
        continue;
    end
    
    for i_ve = 1:size(vessels.cluster_centres,1)
        
        if ~strcmpi(vessels.majority_shapes{i_ve}, 'Distal')
            num_hits = length(unique(vessels.cluster_members{i_ve}));
            
            for i_sh = 1:4
                if strcmpi(vessels.majority_shapes{i_ve}, vessel_shapes{i_sh})
                    vessel_counts2_shape(num_hits,i_sh) = vessel_counts2_shape(num_hits,i_sh) + 1;
                    break;
                end
            end
            for i_sz = 1:5
                if strcmpi(vessels.majority_sizes{i_ve}, vessel_sizes{i_sz})
                    vessel_counts2_size(num_hits,i_sz) = vessel_counts2_size(num_hits,i_sz) + 1;
                    break;
                end
            end
            for i_gr = 1:7
                if strcmpi(vessels.majority_grade, grades{i_gr})
                    vessel_counts2_grade(num_hits,i_gr) = vessel_counts2_grade(num_hits,i_gr) + 1;
                    break;
                end
            end
                
        end
    end
end
display(vessel_counts2_shape);
display(vessel_counts2_size);
display(vessel_counts2_grade);
%%
vessel_counts3 = zeros(3,1);
for i_im = find(markers_per_image == 3)'
    im_num = test_list(i_im).name(1:6);  
    load([cluster_dir im_num '_apex_clusters.mat'], 'vessels');
    if (ismember(grades{6}, vessels.grades) || ismember(grades{7}, vessels.grades))
        continue;
    end
    
    for i_ve = 1:size(vessels.cluster_centres,1)
        
        if ~strcmpi(vessels.majority_shapes{i_ve}, 'Distal')
            num_hits = length(unique(vessels.cluster_members{i_ve}));
            vessel_counts3(num_hits) = vessel_counts3(num_hits) + 1;
        end
    end
end
display(vessel_counts3);

%%
for i_im = find(markers_per_image == 3)'
    im_num = test_list(i_im).name(1:6);  
    load([cluster_dir im_num '_apex_clusters.mat'], 'vessels');
    
    im = u_load([image_dir im_num '.mat']);
    
    figure; imgray(im);
    
    for i_ve = 1:size(vessels.cluster_centres,1)
        
        if ~strcmpi(vessels.majority_shapes{i_ve}, 'Distal')
            plot(vessels.cluster_xy{i_ve}(:,1), vessels.cluster_xy{i_ve}(:,2), 'x');
            plot(vessels.cluster_centres(i_ve,1), vessels.cluster_centres(i_ve,2), '+', 'MarkerSize', 12);
        end
    end
end
%%
total_bad = 0;
for i_im = find(markers_per_image == 3)'
    im_num = test_list(i_im).name(1:6);  
    load([cluster_dir im_num '_apex_clusters.mat'], 'vessels');
    
    total_bad = total_bad + ...
        (ismember(grades{6}, vessels.grades) || ismember(grades{7}, vessels.grades));
end
%%
prob_dir = [nailfoldroot 'data/' rsa_dir 'test/predictions/detection/rf_classification/257273/'];
dist_thresh = 100;
patch_sz2 = 50;

for i_im = 11:30;
    
    if markers_per_image(i_im) > 1
        im_num = test_list(i_im).name(1:6);  
        
        im = u_load([image_dir im_num '.mat']);
        vessel_prob = u_load([prob_dir im_num '_pred.mat']);
        
        load([cluster_dir im_num '_apex_clusters.mat'], 'vessels');

        [vessels] = post_merge_apexes(vessels, vessel_prob, dist_thresh, patch_sz2,...
            0.5, 20, 0);
        
        vessel_idx = sub2ind(size(im), round(vessels.cluster_centres(:,2)), round(vessels.cluster_centres(:,1)));
        
        if ~isempty(vessel_idx)
            fh = figure; imgray(im);
            caxis([min(im(vessel_idx))-10, max(im(vessel_idx))+10])
            plot_apex_clusters(vessels, fh);
        end
    end
end
%%
