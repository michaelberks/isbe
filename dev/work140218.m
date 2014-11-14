%%
num_images = length(image_id_data.im_names);
num_markers = length(image_id_data.marker_list);
marker_counts = zeros(num_markers, 2);
for i_im = 1:num_images
    marked_by = ismember(1:num_markers, image_id_data.marker_idx{i_im});
    if image_id_data.dataset_idx.validation(i_im)
        marker_counts(marked_by,1) = marker_counts(marked_by,1) + 1;
    elseif image_id_data.dataset_idx.test(i_im)
        marker_counts(marked_by,2) = marker_counts(marked_by,2) + 1;
    end
end
%%
image_marked = zeros(5, 2);
for i_im = 1:num_images
    num_marks = length(unique(image_id_data.marker_idx{i_im}));
    if image_id_data.dataset_idx.validation(i_im)
        image_marked(num_marks,1) = image_marked(num_marks,1) + 1;
    elseif image_id_data.dataset_idx.test(i_im)
        image_marked(num_marks,2) = image_marked(num_marks,2) + 1;
    end
end
%%
image_marked_np = zeros(5, 2);
for i_im = 1:num_images
    
    u = unique(image_id_data.marker_idx{i_im});
    u(u == 11) = [];
    num_marks = length(u) + 1;
    
    if image_id_data.dataset_idx.validation(i_im)
        image_marked_np(num_marks,1) = image_marked_np(num_marks,1) + 1;
    elseif image_id_data.dataset_idx.test(i_im)
        image_marked_np(num_marks,2) = image_marked_np(num_marks,2) + 1;
    end
end
%%
%Double marked by category...
marked_twice = zeros(4,2);
for i_im = 1:num_images
    num_marks = length(unique(image_id_data.marker_idx{i_im}));
    
    if num_marks > 1
        if image_id_data.dataset_idx.validation(i_im)
            switch image_id_data.category{i_im}
                case 'S'
                    marked_twice(1,1) = marked_twice(1,1) + 1;
                case 'HC'
                    marked_twice(2,1) = marked_twice(2,1) + 1;
                case 'P'
                    marked_twice(3,1) = marked_twice(3,1) + 1;
                case 'U'
                    marked_twice(4,1) = marked_twice(4,1) + 1;
            end
            
        elseif image_id_data.dataset_idx.test(i_im)
            switch image_id_data.category{i_im}
                case 'S'
                    marked_twice(1,2) = marked_twice(1,2) + 1;
                case 'HC'
                    marked_twice(2,2) = marked_twice(2,2) + 1;
                case 'P'
                    marked_twice(3,2) = marked_twice(3,2) + 1;
                case 'U'
                    marked_twice(4,2) = marked_twice(4,2) + 1;
            end
        end
    end
end
%%

im_repeats = [0 0];
marker_repeats = zeros(num_markers, 2);

for i_im = 1:num_images
    m = image_id_data.marker_idx{i_im};
    [uni_m uni_i] = unique(m);
    n = length(m);
    if length(uni_m) < n
        repeated_by = m(setdiff(1:n, uni_i));
        
        if length(unique(repeated_by)) < length(repeated_by)
            display(['Image ' num2str(i_im) ' repreated more than once: ' num2str(image_id_data.marker_idx{i_im})]);
        end
        
        if image_id_data.dataset_idx.validation(i_im)
            im_repeats(1) = im_repeats(1) + 1;
            marker_repeats(repeated_by,1) = marker_repeats(repeated_by,1) + 1;
        elseif image_id_data.dataset_idx.test(i_im)
            im_repeats(2) = im_repeats(2) + 1;
            marker_repeats(repeated_by,2) = marker_repeats(repeated_by,2) + 1;
        end
    end
end

%%
cluster_repeat_markings(...
    'selected_images',      image_id_data.dataset_idx.validation,...
    'image_id_data',        image_id_data,...
    'test_dir',             'test_half',...
    'data_dir',             [nailfoldroot 'data/rsa_study/'],...
    'cluster_dir',          'repeat_apex_clusters',...
    'marker_list',          [],...
    'discard_repeats',      1,...
    'min_apex_dist',        20,... %The minimum distance allowed between separate clusters (the larger out of this value or the apex width is used for each cluster)
    'include_nondistal',    1,... %Record the position of non-distal vessels, and merge these with vessel marked as distal by a different observer if necessary
    'vessel_prob_dir',      'predictions/detection/rf_classification/296655/',...
    'dist_thresh',          100,... %What's the largest separation allowed between cluster centres
    'patch_sz2',            50,... %What's the halfwidth of the patch in which we look for connectivity
    'connect_thresh',       0.5,... %What's the threshold for deciding 2 clusters are connected?
    'n_connect_pts',        20,...%How many points (spaced between 0 and 1) to test the connectivity of pairs);
    'reduction_factor',     2 ... %Should we resize the ground truth (because we've resize images)
); 
%%
cluster_repeat_markings(...
    'selected_images',      image_id_data.dataset_idx.test,...
    'image_id_data',        image_id_data,...
    'test_dir',             'final_test',...
    'data_dir',             [nailfoldroot 'data/rsa_study/'],...
    'cluster_dir',          'repeat_apex_clusters',...
    'marker_list',          [],...
    'discard_repeats',      1,...
    'min_apex_dist',        20,... %The minimum distance allowed between separate clusters (the larger out of this value or the apex width is used for each cluster)
    'include_nondistal',    1,... %Record the position of non-distal vessels, and merge these with vessel marked as distal by a different observer if necessary
    'vessel_prob_dir',      'predictions/detection/rf_classification/296655/',...
    'dist_thresh',          100,... %What's the largest separation allowed between cluster centres
    'patch_sz2',            50,... %What's the halfwidth of the patch in which we look for connectivity
    'connect_thresh',       0.5,... %What's the threshold for deciding 2 clusters are connected?
    'n_connect_pts',        20,...%How many points (spaced between 0 and 1) to test the connectivity of pairs);
    'reduction_factor',     2 ... %Should we resize the ground truth (because we've resize images)
); 
%%
v_list = dir('C:\isbe\nailfold\data\rsa_study\final_test\repeat_apex_clusters\tmoore\*.mat');
for i_im = 1:20
    im_name = v_list(i_im).name(1:6);
    
    nailfold = u_load(['C:\isbe\nailfold\data\rsa_study\final_test\images\' im_name '.mat']);
    load(['C:\isbe\nailfold\data\rsa_study\final_test\repeat_apex_clusters\tmoore\' v_list(i_im).name]);
    figure; imgray(nailfold);
    plot_apex_clusters(vessels, gcf, 2);
end
%%
apex_marked_once = 0;
apex_marked_twice = 0;
for i_im = 1:20
    load(['C:\isbe\nailfold\data\rsa_study\final_test\repeat_apex_clusters\tmoore\' v_list(i_im).name]);
    for i_ve = 1:length(vessels.cluster_members)
        if ~all(strcmpi('NonDistal', vessels.cluster_shapes{i_ve}))
            if length(vessels.cluster_shapes{i_ve}) > 1
                apex_marked_twice = apex_marked_twice + 1;
            else
                apex_marked_once = apex_marked_once + 1;
            end
        end
    end
end
%%
no_marks  = 0;
for i_im = 1:num_images
    if isempty(image_id_data.markers{i_im})
        no_marks = no_marks + 1;
    end
end