project_dir = 'P:\isbe\nailfold\tom_project\';
root_dir = 'C:\isbe\nailfold\data\rsa_study\';
load('C:\isbe\nailfold\data\rsa_study\data_lists\image_id_data.mat');

im_list = dir([project_dir 'markup\tsaddington\*.txt']);
n_ims = length(im_list);
im_names = cell(n_ims,1);
selected_ims = false(length(image_id_data.im_names),1);

for i_im = 1:n_ims
    im_names{i_im} = im_list(i_im).name(1:6);
    selected_ims = selected_ims | strcmpi(image_id_data.im_names, im_names{i_im});
    
    %[markup] = read_markup_from([project_dir 'markup\tsaddington\' im_list(i_im).name]);
    %write_markup_to(markup, [project_dir 'markup\tsaddington\' im_list(i_im).name], 0.5);
    
end

tom_im_data.im_names = image_id_data.im_names(selected_ims);
tom_im_data.markers = image_id_data.markers(selected_ims,:);
tom_im_data.marker_files = image_id_data.marker_files(selected_ims,:);
tom_im_data.marker_idx = image_id_data.marker_idx(selected_ims,:);
tom_im_data.category = image_id_data.category(selected_ims,:);
%%
for i_im = 1:n_ims
    tom_im_data.markers{i_im}(end+1) = {'tsaddington'};
    m_name = dir([project_dir 'markup\tsaddington\*' tom_im_data.im_names{i_im} '*.txt']);
    tom_im_data.marker_files{i_im}(end+1) = { [project_dir 'markup\tsaddington\' m_name(1).name]};
    tom_im_data.marker_idx{i_im}(end+1) = 14;
end
%%
markers = {...
    'aherrick'
    'asulli'
    'fravera'
    'hhofstee'
    'jallen'
    'khowell'
    'manderson'
    'mcutolo'
    'mwildt'
    'ppyrkotsch'
    'rhesselstrand'
    'tmoore'
    'vsmith'
    'tsaddington'};

for i_ma = 1:length(markers)-1
    create_folder([project_dir 'markup\' markers{i_ma}]);
end
create_folder([project_dir 'predictions\detection\rf_classification\296655\']);
create_folder([project_dir 'images\']);


for i_im = 1:n_ims
    
    if ~exist([project_dir 'images\' im_names{i_im} '.mat'], 'file') &&...
        exist([root_dir 'master_set\images\' im_names{i_im} '.mat'], 'file')
        copyfile(...
            [root_dir 'master_set\images\' im_names{i_im} '.mat'],...
            [project_dir 'images\' im_names{i_im} '.mat']);
    end
    
    if ~exist([project_dir 'predictions\detection\rf_classification\296655\' im_names{i_im} '_pred.mat'], 'file') &&...
        exist([root_dir 'master_set\predictions\detection\rf_classification\296655\' im_names{i_im} '_pred.mat'], 'file')
        copyfile(...
            [root_dir 'master_set\predictions\detection\rf_classification\296655\' im_names{i_im} '_pred.mat'],...
            [project_dir 'predictions\detection\rf_classification\296655\' im_names{i_im} '_pred.mat']);
    end
    
    for i_ma = 1:length(markers)-1
    
        m_list = dir([root_dir 'markup\' markers{i_ma} '\*' im_names{i_im} '*.txt']);
        
        for i_mf = 1:length(m_list)
            copyfile(...
                [root_dir 'markup\' markers{i_ma} '\' m_list(i_mf).name],...
                [project_dir 'markup\' markers{i_ma} '\' m_list(i_mf).name]);
        end
    end
end
%%
cluster_vessel_apices_set(...
    'selected_images',      true(30,1), ...
    'image_id_data',        tom_im_data,...
    'test_dir',             [],...
    'data_dir',             project_dir,...
    'cluster_dir',          'apex_clusters',...
    'post_cluster_dir',     [],...
    'apex_gt_dir',          'apex_gt',...
    'marker_list',          markers,...
    'min_apex_dist',        20,... %The minimum distance allowed between separate clusters (the larger out of this value or the apex width is used for each cluster)
    'include_nondistal',    1,... %Record the position of non-distal vessels, and merge these with vessel marked as distal by a different observer if necessary
    'vessel_prob_dir',      'predictions/detection/rf_classification/296655/',...
    'dist_thresh',          100,... %What's the largest separation allowed between cluster centres
    'patch_sz2',            50,... %What's the halfwidth of the patch in which we look for connectivity
    'connect_thresh',       0.5,... %What's the threshold for deciding 2 clusters are connected?
    'n_connect_pts',        20,...%How many points (spaced between 0 and 1) to test the connectivity of pairs);
    'reduction_factor',     2, ... %Should we resize the ground truth (because we've resize images)
    'do_post_merge',        0 ...
); 
%%
ap_list = dir(['P:\isbe\nailfold\tom_project\apex_clusters\*.mat']);
for i_ap = 1:30
    nailfold = imread([project_dir 'images\' im_names{i_ap} '.png']);
    load(['P:\isbe\nailfold\tom_project\apex_clusters\' ap_list(i_ap).name]);
    f = figure;
    imgray(nailfold);
    plot_apex_clusters(vessels, f);
end 
%%
full_exists = false(30,1);
for i_im = 1:30
    full_exists(i_im) = exist([root_dir 'images\' im_names{i_im} '.png'], 'file');
    copyfile(...
        [root_dir 'images\' im_names{i_im} '.png'],...
        [project_dir 'images\' im_names{i_im} '.png']);
end
%%
markers = {...
    'aherrick'
    'asulli'
    'fravera'
    'hhofstee'
    'jallen'
    'khowell'
    'manderson'
    'mcutolo'
    'mwildt'
    'ppyrkotsch'
    'rhesselstrand'
    'tmoore'
    'vsmith'
    'tsaddington'};

for i_ma = 1:length(markers)
    mf_list = dir([project_dir 'markup\' markers{i_ma} '\*.txt']);
    for i_mf = 1:length(mf_list)
        vessel_data_file = [project_dir 'markup\' markers{i_ma} '\' mf_list(i_mf).name];
        excel_file = [project_dir 'markup\' markers{i_ma} '\' mf_list(i_mf).name(1:end-3) 'xls'];
        read_vessels_to_excel(vessel_data_file, excel_file, 1);
    end
end
%%
ap_list = dir('C:\isbe\nailfold\data\rsa_study\master_set\apex_clusters\*.mat');
grade_counts = zeros(7,4);
for i_ap = 1:length(ap_list)
    load(['C:\isbe\nailfold\data\rsa_study\master_set\apex_clusters\' ap_list(i_ap).name]);
    m_idx = find(vessels.markers == 1, 1);
    if ~isempty(m_idx)
        grade = vessels.grades{m_idx};
        im_name = ap_list(i_ap).name(1:6);
        im_idx = strcmpi(image_id_data.im_names, im_name);
        category = image_id_data.category{im_idx};
        
        switch category
            
            case 'HC'
                c_idx = 1;
            case 'P'
                c_idx = 2;
            case 'S'
                c_idx = 3;
            case 'U'
                c_idx = 4;
                
            otherwise
                continue;
        end
        
        
        switch grade
            case 'Normal'
                g_idx = 1;
            case 'Early'
                g_idx = 2;
            case 'Active'
                g_idx = 3;
            case 'Late'
                g_idx = 4;
            case 'Non-specific'
                g_idx = 5;
            case 'Ungradeable_Quality'
                g_idx = 6;
            case 'Ungradeable_Condition'
                g_idx = 7;
                
            otherwise
                continue;
        end
        
        grade_counts(g_idx, c_idx) = grade_counts(g_idx, c_idx) + 1;
    end
end

