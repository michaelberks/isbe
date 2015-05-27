load('C:\isbe\nailfold\data\rsa_study\data_lists\image_id_data.mat');
%%
cluster_vessel_apices_set(...
    'selected_images',      ~(image_id_data.dataset_idx.validation | image_id_data.dataset_idx.test) & (1:length(image_id_data.dataset_idx.validation)>=1936)', ...
    'image_id_data',        image_id_data,...
    'test_dir',             'master_set',...
    'data_dir',             [nailfoldroot 'data/rsa_study/'],...
    'cluster_dir',          'apex_clusters_with_repeats',...
    'post_cluster_dir',     'apex_clusters_with_repeats_merged',...
    'apex_gt_dir',          'apex_gt',...
    'marker_list',          [],...
    'discard_repeats',      0,...
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
cluster_dir = 'C:\isbe\nailfold\data\rsa_study\master_set\apex_clusters_with_repeats_merged\';
im_list = dir([cluster_dir '*.mat']);

all_vessels_long_form = cell(1e4,11);
curr_vess = 1;
for i_im = 1:length(im_list);
    vessels = u_load([cluster_dir im_list(i_im).name]);
    
    im_name = im_list(i_im).name(1:6);
    
    num_vessels = length(vessels.cluster_members);
    if num_vessels == 0;
        display(['Skipping ' im_name ', no vessels marked']);
        continue;
    end
    
    for i_ve = 1:num_vessels
        
        for i_ma = 1:length(vessels.cluster_members{i_ve})
        
            all_vessels_long_form{curr_vess,1} = i_im; %Image #
            all_vessels_long_form{curr_vess,2} = im_name; %Image name
            all_vessels_long_form{curr_vess,3} = i_ve; %Vessel #

            observer = vessels.cluster_members{i_ve}(i_ma);
            if observer < 100
                all_vessels_long_form{curr_vess,4} = observer; %Observer #
                all_vessels_long_form{curr_vess,5} = 1; %Repeat #
            elseif vessels.cluster_members{i_ve}(i_ma) < 200
                all_vessels_long_form{curr_vess,4} = observer-100; %Observer #
                all_vessels_long_form{curr_vess,5} = 2; %Repeat #
            elseif vessels.cluster_members{i_ve}(i_ma) < 300
                all_vessels_long_form{curr_vess,4} = observer-200; %Observer #
                all_vessels_long_form{curr_vess,5} = 3; %Repeat #
            elseif vessels.cluster_members{i_ve}(i_ma) < 400
                all_vessels_long_form{curr_vess,4} = observer-300; %Observer #
                all_vessels_long_form{curr_vess,5} = 4; %Repeat #
            else
                all_vessels_long_form{curr_vess,4} = observer-400; %Observer #
                all_vessels_long_form{curr_vess,5} = 5; %Repeat #
            end
            
            if ~strcmp(vessels.cluster_shapes{i_ve}{i_ma}, 'NonDistal')
                all_vessels_long_form{curr_vess,6} = 1;
                all_vessels_long_form{curr_vess,7} = vessels.cluster_widths{i_ve}(i_ma);
            else
                all_vessels_long_form{curr_vess,6} = 0;
                all_vessels_long_form{curr_vess,7} = -1;
            end
            
            %Shape and size descriptors
            all_vessels_long_form{curr_vess,8} = vessels.cluster_shapes{i_ve}{i_ma};
            all_vessels_long_form{curr_vess,9} = vessels.cluster_shapes{i_ve}{i_ma};
            
            %Anchor point
            all_vessels_long_form{curr_vess,10} = vessels.cluster_xy{i_ve}(i_ma,1); %Anchor x
            all_vessels_long_form{curr_vess,11} = vessels.cluster_xy{i_ve}(i_ma,2);%Achor y
            
            curr_vess = curr_vess + 1;
        end
    end
end

%Discard any additional rows from the pre-allocated container
all_vessels_long_form(curr_vess:end,:) = [];
        
    
    
    