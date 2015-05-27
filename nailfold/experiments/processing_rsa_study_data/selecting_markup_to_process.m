%--------------------------------------------------------------------------
% Selecting vessel training data, this is fully annotated, with vessel
% centre paths and edges marked and is used to train the pixel-wise vessel
% detection, orientation and width RF predictors
%--------------------------------------------------------------------------

%All the images and markup files are original stored on a Faculty server
%that can be accessed using the details below...
%Host: sftp.rss.mhs.man.ac.uk
%port 22
%username: ptran
%protocol: sftp
%allow fallback? Yes
%p:nailfold2000


%--------------------------------------------------------------------------
%Four stage sorting of vessel markups
%--------------------------------------------------------------------------
% 1) Quality score each markup as reject, skip or process
% 2) Quality score all the process vessels as good or not
% 3) Manually markup inner and outer edges of the good vessels
% 4) Mark a field-of-view mask for each annotated good vessel
%--------------------------------------------------------------------------

%Repeat for each class of vessel size: normal, enlarged, giant
vessel_size = 'giant';
vessel_dir = ['C:\isbe\nailfold\data\rsa_study\apexes\' vessel_size '\'];
v_list = dir([vessel_dir, '*vessel.mat']);

%% 1) Quality score each markup as reject, skip or process
quick_vessel_segmentation(vessel_dir);
%%

skip_count = 0;
reject_count = 0;
process_count = 0;

process_selection = false(1, length(v_list));
for i_ve = 1:length(v_list)
    load([vessel_dir v_list(i_ve).name]);
    
    switch vessel_struc.quality
        case 0
            reject_count = reject_count + 1;
        case 1
            skip_count = skip_count + 1;
        case 2
            process_count = process_count + 1;
            process_selection(i_ve) = 1;
    end
end
process_selection = find(process_selection);

%% 2) Which of the process vessels are really good ? i.e. can have both edges annotated
quick_vessel_segmentation2(vessel_dir, process_selection);
%%
good_selection = false(1, length(v_list));
for i_ve = 1:length(v_list)
    load([vessel_dir v_list(i_ve).name]);
    
    switch vessel_struc.quality
        case 3
            good_selection(i_ve) = 1;
    end
end
good_selection = find(good_selection);
%% 3) Mark up both edges of the good selection
contour_dir = ['C:\isbe\nailfold\data\rsa_study\vessel_contours\' vessel_size '\'];
create_folder(contour_dir);
quick_vessel_segmentation3(vessel_dir, contour_dir, good_selection);

%Convert each edge markup into a proper contour (with matched edges and a
%tue centre line)
convert_edge_markup_batch(contour_dir, [], 1);
%% 4) Mark up a fov mask for each good vessel
fov_dir = ['C:\isbe\nailfold\data\rsa_study\fov_masks\' vessel_size '\'];
create_folder(fov_dir);
quick_fov_markup(vessel_dir, fov_dir, good_selection);

%% 5) Mark any apexes on the centreline
quick_apex_marking(contour_dir);

%% 6) Compute some statistics on the widths of the vessels and work out which images they were extracted from
vessel_sizes = {'normal', 'enlarged', 'giant'};
vessel_lengths = cell(1,3);
vessel_apex_widths_orig = cell(1,3);
vessel_apex_widths_new = cell(1,3);
vessel_apex_shapes = cell(1,3);
vessel_names = cell(1,3);
all_vessel_widths = cell(1,3);

for i_sz = 1:3
    vessel_dir = ['C:\isbe\nailfold\data\rsa_study\apexes\' vessel_sizes{i_sz} '\'];
    contour_dir = ['C:\isbe\nailfold\data\rsa_study\vessel_contours\' vessel_sizes{i_sz} '\'];
    v_files = dir([contour_dir '*contour.mat']);
    
    num_vessels = length(v_files);
    
    
    vessel_lengths{i_sz} = zeros(num_vessels,1); %a
    vessel_apex_widths_orig{i_sz} = zeros(num_vessels,1); %b
    vessel_apex_widths_new{i_sz} = zeros(num_vessels,1); %c
    vessel_apex_shapes{i_sz} = cell(num_vessels,1); %d
    vessel_names{i_sz} = cell(num_vessels,1); %e   
    all_vessel_widths{i_sz} = []; %f
    
    for i_v = 1:length(v_files)
        %load data
        contour_struc = load([contour_dir v_files(i_v).name]);
        apex_struc = load([vessel_dir v_files(i_v).name(1:8) '.mat']);
        vessel_struc = u_load([vessel_dir v_files(i_v).name(1:8) '_vessel.mat']);
            
        if ~isfield(vessel_struc, 'vessel_properties')
            continue;
        end
        
        [~, im_name] = fileparts(apex_struc.apex_properties.nailfold_name);
        apex_xy(:,1) = apex_struc.apex_xy(:,1) +...
            apex_struc.apex_properties.sc - vessel_struc.vessel_properties.sc;
        apex_xy(:,2) = apex_struc.apex_xy(:,2) +...
            apex_struc.apex_properties.sr - vessel_struc.vessel_properties.sr;

        %Compute centre and widths from contour edges
        vessel_centre = (contour_struc.inner_edge + contour_struc.outer_edge)/2;
        v_widths_i = sqrt(sum((contour_struc.inner_edge - contour_struc.outer_edge).^2,2));
        
        %Workout which centre point lies closest to the originally marked
        %apex
        dists = sum(bsxfun(@minus, vessel_centre, mean(apex_xy)).^2,2);
        [~, min_i] = min(dists);
        
        %Set a
        vessel_lengths{i_sz}(i_v) = sum(sqrt(sum(diff(vessel_centre).^2,2)));
        %Set b
        vessel_apex_widths_orig{i_sz}(i_v) = sqrt(sum(diff(apex_xy).^2,2));
        %Set c
        vessel_apex_widths_new{i_sz}(i_v) = v_widths_i(min_i);
        %Set d
        vessel_apex_shapes{i_sz}{i_v} = apex_struc.apex_properties.vessel_shape;
        %Set e
        vessel_names{i_sz}{i_v} = im_name;
        %Set f        
        all_vessel_widths{i_sz} = [all_vessel_widths{i_sz}; v_widths_i];   
        
        if 0%vessel_apex_widths_orig{i_sz}(i_v) > 64
            if isfield(vessel_struc, 'vessel_patch');
                figure; imgray(vessel_struc.vessel_patch);
                title(vessel_sizes{i_sz});
            end
        end
            
    end
end

counts(1,:) = hist(all_vessel_widths{1}, 1:2:250);
counts(2,:) = hist(all_vessel_widths{2}, 1:2:250);
counts(3,:) = hist(all_vessel_widths{3}, 1:2:250);
figure; plot(1:2:250, bsxfun(@rdivide, counts, sum(counts,2))', 'linewidth', 2);

save('C:\isbe\nailfold\data\rsa_study\training\training_data_stats.mat',...
    'vessel_lengths', 'vessel_apex_widths_orig', 'vessel_apex_widths_new', 'vessel_apex_shapes', 'vessel_names', 'all_vessel_widths');
%% ------------------------------------------------------------------------
%% 
% Now compute new vessel apex points given the selected, annotated vessels
vessel_sizes = {'normal', 'enlarged', 'giant'};
v_step_sz = 2;
num_v_pts = 15;
v_sample_pts = v_step_sz*(-num_v_pts:num_v_pts);

mean_apex_width =  mean([vessel_apex_widths_new{1};vessel_apex_widths_new{2}]);

for i_sz = 1:2
    vessel_dir = ['C:\isbe\nailfold\data\rsa_study\apexes\' vessel_sizes{i_sz} '\'];
    contour_dir = ['C:\isbe\nailfold\data\rsa_study\vessel_contours\' vessel_sizes{i_sz} '\'];
    v_files = dir([contour_dir '*contour.mat']);
    
    num_vessels = length(v_files);
       
    for i_v = 1:length(v_files)
        %load data
        contour_struc = load([contour_dir v_files(i_v).name]);
        apex_struc = load([vessel_dir v_files(i_v).name(1:8) '.mat']);
        vessel_struc = u_load([vessel_dir v_files(i_v).name(1:8) '_vessel.mat']);
        
        apex_xy(:,1) = apex_struc.apex_xy(:,1) +...
            apex_struc.apex_properties.sc - vessel_struc.vessel_properties.sc;
        apex_xy(:,2) = apex_struc.apex_xy(:,2) +...
            apex_struc.apex_properties.sr - vessel_struc.vessel_properties.sr;
        
        %Compute centre and widths from contour edges
        vessel_centre = (contour_struc.inner_edge + contour_struc.outer_edge)/2;
        v_widths_i = sqrt(sum((contour_struc.inner_edge - contour_struc.outer_edge).^2,2));
         
        %Workout which centre point lies closest to the originally marked
        %apex and use this to work out the apex width
        dists = sum(bsxfun(@minus, vessel_centre, mean(apex_xy)).^2,2);
        [~, apex_i] = min(dists);
        apex_width = v_widths_i(apex_i);
        
        %Now select the vessel top
        
        scale_factor = apex_width / mean_apex_width;
        v_sample_pts_i = v_sample_pts * scale_factor;
        
        v_running_dist = [0; cumsum(sqrt(sum(diff(vessel_centre).^2,2)))];
        v_sample_dists = v_running_dist(apex_i) + v_sample_pts_i;
        v_pts = interp1(v_running_dist, vessel_centre, v_sample_dists, 'linear');
        
        vessel_struc.v_pts_new = v_pts;
        save([vessel_dir v_files(i_v).name(1:8) '_vessel.mat'], 'vessel_struc');
        
        if i_v <= 5
            figure; imgray(vessel_struc.vessel_patch);
            plot(v_pts(:,1), v_pts(:,2), 'g');
            plot(v_pts(:,1), v_pts(:,2), 'r.');  
            plot(apex_xy(:,1), apex_xy(:,2), 'rx');
            plot(apex_xy(:,1), apex_xy(:,2), 'y');
        end
    end
end
%% -------------------------------------------------------------------------
%% 7) Now make vessel masks and copy the various images and maps to specifc sets of training data
%
%Repeat for each class of vessel size: normal, enlarged, giant
vessel_sizes = {'normal', 'enlarged', 'giant'};
set_dirs = {'set12g'};

set_split = 1; %set to 0.5 to split into a separate training/test set
[uni_names uni_idx] = unique([vessel_names{1}; vessel_names{2}]);
set_names = {uni_names(1:floor(end*set_split)), uni_names(ceil(end*set_split):end)};

for i_st = 1
    vessel_mask_dir = ['C:\isbe\nailfold\data\rsa_study\' set_dirs{i_st} '\vessel_masks\'];
    vessel_cmask_dir = ['C:\isbe\nailfold\data\rsa_study\' set_dirs{i_st} '\vessel_centre_masks\'];
    vessel_image_dir = ['C:\isbe\nailfold\data\rsa_study\' set_dirs{i_st} '\images\'];
    vessel_imagen_dir = ['C:\isbe\nailfold\data\rsa_study\' set_dirs{i_st} '\images_n\'];
    vessel_width_dir = ['C:\isbe\nailfold\data\rsa_study\' set_dirs{i_st} '\width_maps\'];
    vessel_ori_dir = ['C:\isbe\nailfold\data\rsa_study\' set_dirs{i_st} '\orientations\'];
    vessel_fov_dir = ['C:\isbe\nailfold\data\rsa_study\' set_dirs{i_st} '\fov_masks\'];
    vessel_contour_dir = ['C:\isbe\nailfold\data\rsa_study\' set_dirs{i_st} '\vessel_contours\'];
    
    create_folder(vessel_mask_dir);
    create_folder(vessel_cmask_dir);
    create_folder(vessel_image_dir);
    create_folder(vessel_imagen_dir);
    create_folder(vessel_width_dir);
    create_folder(vessel_ori_dir);
    create_folder(vessel_fov_dir);
    create_folder(vessel_contour_dir);
end
for i_sz = 1:3

    vessel_dir = ['C:\isbe\nailfold\data\rsa_study\apexes\' vessel_sizes{i_sz} '\'];
    contour_dir = ['C:\isbe\nailfold\data\rsa_study\vessel_contours\' vessel_sizes{i_sz} '\'];
    
    make_masks = 1;
    copy_images = 1;
    copy_masks = 1;

    v_files = dir([contour_dir '*contour.mat']);
    for i_v = 1:length(v_files)
        contour_struc = load([contour_dir v_files(i_v).name]);
        vessel_struc = u_load([vessel_dir v_files(i_v).name(1:end-12) '.mat']);

        if ismember(vessel_names{i_sz}{i_v}, set_names{1})
            i_st = 1;
        else
            i_st = 2;
        end
        
        vessel_mask_dir = ['C:\isbe\nailfold\data\rsa_study\' set_dirs{i_st} '\vessel_masks\'];
        vessel_cmask_dir = ['C:\isbe\nailfold\data\rsa_study\' set_dirs{i_st} '\vessel_centre_masks\'];
        vessel_image_dir = ['C:\isbe\nailfold\data\rsa_study\' set_dirs{i_st} '\images\'];
        vessel_imagen_dir = ['C:\isbe\nailfold\data\rsa_study\' set_dirs{i_st} '\images_n\'];
        vessel_width_dir = ['C:\isbe\nailfold\data\rsa_study\' set_dirs{i_st} '\width_maps\'];
        vessel_ori_dir = ['C:\isbe\nailfold\data\rsa_study\' set_dirs{i_st} '\orientations\'];
        vessel_fov_dir = ['C:\isbe\nailfold\data\rsa_study\' set_dirs{i_st} '\fov_masks\'];
        vessel_contour_dir = ['C:\isbe\nailfold\data\rsa_study\' set_dirs{i_st} '\vessel_contours\'];
    
        if make_masks
            [rows cols] = size(vessel_struc.vessel_patch);
            [vessel_mask vessel_centre_mask width_map ori_map] =...
                make_mask_from_contour(contour_struc.outer_edge, contour_struc.inner_edge, rows, cols);
            save([vessel_mask_dir vessel_sizes{i_sz} v_files(i_v).name(1:end-12) '_v_mask.mat'], 'vessel_mask');
            save([vessel_cmask_dir vessel_sizes{i_sz} v_files(i_v).name(1:end-12) '_v_cmask.mat'], 'vessel_centre_mask');
            save([vessel_width_dir vessel_sizes{i_sz} v_files(i_v).name(1:end-12) '_width.mat'], 'width_map');
            save([vessel_ori_dir vessel_sizes{i_sz} v_files(i_v).name(1:end-12) '_ori.mat'], 'ori_map');
            copyfile([contour_dir v_files(i_v).name],...
                [vessel_contour_dir vessel_sizes{i_sz} v_files(i_v).name(1:end-12) '_vc.mat']);

            if i_v <= 10
                figure; 
                subplot(2,3,1); imgray(vessel_centre_mask);
                subplot(2,3,2); imgray(width_map);
                subplot(2,3,3); imgray(complex2rgb(ori_map));
            end

        end
        if copy_images
            vessel_patch = double(vessel_struc.vessel_patch);
            g = gaussian_filters_1d(32, 96);
            g = g / sum(g);
            patch_mask = make_nailfold_mosaic_mask(vessel_patch);
            im_edges = conv2(g', g, double(patch_mask), 'same');
            vessel_patch_smoothed = conv2(g', g, vessel_patch, 'same') ./ im_edges;
            vessel_patch_equalised = vessel_patch - vessel_patch_smoothed;
            vessel_patch_equalised(~patch_mask) = 0;
            if i_v <= 10 
                subplot(2,3,1); imgray(vessel_patch);
                subplot(2,3,2); imgray(vessel_patch_equalised);
                subplot(2,3,3); imgray(patch_mask);
            end
            save([vessel_image_dir vessel_sizes{i_sz} v_files(i_v).name(1:end-12) '.mat'], 'vessel_patch');
            save([vessel_imagen_dir vessel_sizes{i_sz} v_files(i_v).name(1:end-12) '.mat'], 'vessel_patch_equalised');
        end
        if copy_masks
            fov_mask = u_load(['C:\isbe\nailfold\data\rsa_study\fov_masks\' vessel_sizes{i_sz} '\' v_files(i_v).name(1:8) '_vessel_fov_mask.mat']);
            save([vessel_fov_dir vessel_sizes{i_sz} v_files(i_v).name(1:end-12) '_f_mask.mat'], 'fov_mask');
        end
    end
end
%% 7) Make a 90 degree copy of all the images
for i_st = 1:2
    base_dir = ['C:\isbe\nailfold\data\rsa_study\set' num2str(i_st) 'g\'];
    im_list = dir([base_dir 'images\*mat']);
    im_types = {...
        'images/', '.mat';...
        'images_n/', '.mat';...
        'fov_masks/', '_f_mask.mat';...
        'vessel_masks/', '_v_mask.mat';
        'vessel_centre_masks/', '_v_cmask.mat';
        'width_maps/', '_width.mat';
        'orientations/', '_ori.mat'};
    num_types = size(im_types,1);
    for i_im = 1:length(im_list);
        im_name = im_list(i_im).name(1:end-4);
        l_or_r = 1 - 2*(rand>.5);

        for i_type = 1:num_types
            vessel_im = u_load([base_dir im_types{i_type,1} im_name im_types{i_type,2}]);
            vessel_im = rot90(vessel_im, l_or_r);

            if strcmpi(im_types{i_type,1}(1), 'o')
                vessel_im = -vessel_im;
                display(['Image ' num2str(i_im) ' flipped orientation']);
            end
            save([base_dir im_types{i_type,1} '90_' im_name im_types{i_type,2}], 'vessel_im');
        end
    end
end
%% 8) Check how many vessel and centre points we have in the training data
cmask_list = dir([vessel_cmask_dir '*.mat']);
vmask_list = dir([vessel_mask_dir '*.mat']);

v_count = 0;
c_count = 0;
for i_ve = 1:length(cmask_list)
    c_mask = u_load([vessel_cmask_dir cmask_list(i_ve).name]);
    v_mask = u_load([vessel_mask_dir vmask_list(i_ve).name]);
    
    c_count = c_count + sum(c_mask(:));
    v_count = v_count + sum(v_mask(:));
end
%% ------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%% ------------------------------------------------------------------------
%% Next we select the the remaining images to go into a test set
% This will actually end being used as a validation set - it was all the
% images we had circa Jan 2013. We have subsequently updated the images,
% and the new set will be used as a genuine test set
%% 1) Copy over all the main nailfold images not annotated into a test directory
normal_names = unique(vessel_names{1});
enlarged_names = unique(vessel_names{2});
all_names = unique([vessel_names{1}; vessel_names{2}]);

main_dir = 'C:\isbe\nailfold\data\rsa_study\';
test_dir = 'C:\isbe\nailfold\data\rsa_study\test\';
mkdir([test_dir 'images']);
mkdir([test_dir 'fov_masks']);

all_images = dir([main_dir 'images\*.png']);

for i_im = 1:length(all_images);
    [~, im_name] = fileparts(all_images(i_im).name);
    if ~ismember(im_name, all_names)
        nailfold = double(imread([main_dir 'images\' all_images(i_im).name]));
        fov_mask = make_nailfold_mosaic_mask(nailfold);
        
        save([test_dir 'images\' im_name '.mat'], 'nailfold');
        save([test_dir 'fov_masks\' im_name '_f_mask.mat'], 'fov_mask');
    end
end

%% 2) Get apex observer ground truth from amrkup files
load('C:\isbe\nailfold\data\rsa_study\data_lists\image_id_data.mat');

cluster_vessel_apices_set(...
    'selected_images',      image_id_data.dataset_idx.validation | image_id_data.dataset_idx.test, ...
    'image_id_data',        image_id_data,...
    'test_dir',             'master_set',...
    'data_dir',             [nailfoldroot 'data/rsa_study/'],...
    'cluster_dir',          'apex_clusters',...
    'post_cluster_dir',     'apex_clusters_merged',...
    'apex_gt_dir',          'apex_gt',...
    'marker_list',          [],...
    'min_apex_dist',        20,... %The minimum distance allowed between separate clusters (the larger out of this value or the apex width is used for each cluster)
    'include_nondistal',    1,... %Record the position of non-distal vessels, and merge these with vessel marked as distal by a different observer if necessary
    'vessel_prob_dir',      'predictions/detection/rf_classification/296655/',...
    'dist_thresh',          100,... %What's the largest separation allowed between cluster centres
    'patch_sz2',            50,... %What's the halfwidth of the patch in which we look for connectivity
    'connect_thresh',       0.5,... %What's the threshold for deciding 2 clusters are connected?
    'n_connect_pts',        20,...%How many points (spaced between 0 and 1) to test the connectivity of pairs);
    'reduction_factor',     2 ... %Should we resize the ground truth (because we've resize images)
); 

%% ------------------------------------------------------------------------
%Everything we get from when the initial test above was fixed we add to
% a new genuine test set
%--------------------------------------------------------------------------
%% 1) Copy, resize the images and make field-of-view masks
load('C:\isbe\nailfold\data\rsa_study\data_lists\image_id_data.mat');
image_dir = 'C:\isbe\nailfold\data\rsa_study\images\';
test_dir = 'C:\isbe\nailfold\data\rsa_study\master_set\';
reduction_factor = 2;

test_names = image_id_data.im_names(image_id_data.dataset_idx.test);

create_folder([test_dir 'images']);
create_folder([test_dir 'fov_masks']);
%
num_images = length(test_names);
for i_im = 1:num_images
    display(['Processing image ' num2str(i_im) ' of ' num2str(num_images)]);
    
    if ~exist([test_dir 'images\' test_names{i_im} '.mat'], 'file')
        nailfold = double(imread([image_dir test_names{i_im} '.png']));
        nailfold = imresize(nailfold, 0.5);

        fov_mask = make_nailfold_mosaic_mask(nailfold, 250, 5);

        save([test_dir 'images\' test_names{i_im} '.mat'], 'nailfold');
        save([test_dir 'fov_masks\' test_names{i_im} '_f_mask.mat'], 'fov_mask');
    end   
end

%% 2) Get apex observer ground truth from markup files
load('C:\isbe\nailfold\data\rsa_study\data_lists\image_id_data.mat');
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
    'vsmith'};

cluster_vessel_apices_set(...
    'selected_images',      image_id_data.dataset_idx.test, ...
    'image_id_data',        image_id_data,...
    'test_dir',             'master_set',...
    'data_dir',             [nailfoldroot 'data/rsa_study/'],...
    'cluster_dir',          'apex_clusters',...
    'post_cluster_dir',     'apex_clusters_merged',...
    'apex_gt_dir',          'apex_gt',...
    'marker_list',          [],...
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
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%% 3) Copying over the final batch of images
image_dir = 'C:\isbe\nailfold\data\rsa_study\images\';
all_caps = dir([image_dir 'all_data\*c.bmp']);
for i_im = 1:length(all_caps)
    [~, im_name] = fileparts(all_caps(i_im).name);
    
    if ~exist([image_dir im_name '.png'], 'file')
        nailfold = imread([image_dir 'all_data\' all_caps(i_im).name]);
        nailfold = nailfold(:,:,1);
        imwrite(nailfold, [image_dir im_name '.png']);
        display(['Copying ' im_name]);        
    end
    delete([image_dir 'all_data\' all_caps(i_im).name]);
end
%%
im_list = dir([base_dir 'images_png\*.png']);
fname = 'C:\isbe\nailfold\data\rsa_study\data_lists\set12g_half_image_list.txt';
fid1 = fopen(fname, 'wt');
for i_im = 1:length(im_list)
    fprintf(fid1,'%s \n', im_list(i_im).name);
end
fclose(fid1);
    