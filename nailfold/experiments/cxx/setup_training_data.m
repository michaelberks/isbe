%Making data for C++

base_dir = 'C:\isbe\nailfold\data\rsa_study\';
matlab_dir = [base_dir 'set12g\'];
cxx_dir = [base_dir 'cxx\'];
vessel_mask_dir = [matlab_dir 'vessel_masks\'];
vessel_cmask_dir = [matlab_dir 'vessel_centre_masks\'];
vessel_image_dir = [matlab_dir 'images\'];
vessel_width_dir = [matlab_dir 'width_maps\'];
vessel_ori_dir = [matlab_dir 'orientations\'];
vessel_fov_dir = [matlab_dir 'fov_masks\'];

vessel_mask_dirc = [cxx_dir 'vessel_masks\'];
vessel_cmask_dirc = [cxx_dir 'vessel_centre_masks\'];
vessel_image_dirc = [cxx_dir 'images\'];
vessel_width_dirc = [cxx_dir 'width_maps\'];
vessel_ori_dirc = [cxx_dir 'orientations\'];
vessel_fov_dirc = [cxx_dir 'fov_masks\'];

vessel_cont_dirc = [cxx_dir 'vessel_contours\'];
    
vessel_sizes = {'normal', 'enlarged', 'giant'};

im_number = 1;
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
        
        full_image = imread(vessel_struc.vessel_properties.nailfold_name);

        [rows cols] = size(full_image);
        [patch_rows patch_cols] = size(vessel_struc.vessel_patch);
        
        sr0 = vessel_struc.vessel_properties.sr;
        er0 = vessel_struc.vessel_properties.sr + patch_rows - 1;
        sc0 = vessel_struc.vessel_properties.sc;
        ec0 = vessel_struc.vessel_properties.sc + patch_cols - 1;
        
        extra_rows = floor((512-patch_rows)/2);
        extra_cols = floor((512-patch_cols)/2);
        if extra_rows <= 0
            sr1 = sr0;
            er1 = er0;
            
        elseif sr0 < extra_rows;
            sr1 = 1;
            er1 = 512;
            
        elseif er0+extra_rows > rows
            er1 = rows;
            sr1 = rows-511;
            
        else
            sr1 = sr0-extra_rows;
            er1 = sr1 + 511;
        end
            
        if extra_cols <= 0
            sc1 = sc0;
            ec1 = ec0;
            
        elseif sc0 < extra_cols;
            sc1 = 1;
            ec1 = 512;
            
        elseif ec0+extra_cols > cols
            ec1 = cols;
            sc1 = cols-511;
            
        else
            sc1 = sc0-extra_cols;
            ec1 = sc1 + 511;
        end
        
        pad_t = sr0 - sr1;
        pad_b = er1 - er0;
        pad_l = sc0 - sc1;
        pad_r = ec1 - ec0;

        vessel_patch = full_image(sr1:er1, sc1:ec1);
        
        %Now pad everything appropriately
        vessel_mask = u_load([vessel_mask_dir vessel_sizes{i_sz} v_files(i_v).name(1:end-12) '_v_mask.mat']);
        vessel_mask = padarray(padarray(vessel_mask, [pad_t pad_l], 0, 'pre'), [pad_b pad_r], 0, 'post');
        
        vessel_centre_mask = u_load([vessel_cmask_dir vessel_sizes{i_sz} v_files(i_v).name(1:end-12) '_v_cmask.mat']);
        vessel_centre_mask = padarray(padarray(vessel_centre_mask, [pad_t pad_l], 0, 'pre'), [pad_b pad_r], 0, 'post');
        
        width_map = u_load([vessel_width_dir vessel_sizes{i_sz} v_files(i_v).name(1:end-12) '_width.mat']);
        width_map = padarray(padarray(width_map, [pad_t pad_l], 0, 'pre'), [pad_b pad_r], 0, 'post');
        
        ori_map = u_load([vessel_ori_dir vessel_sizes{i_sz} v_files(i_v).name(1:end-12) '_ori.mat']);
        ori_map = padarray(padarray(ori_map, [pad_t pad_l], 0, 'pre'), [pad_b pad_r], complex(0,0), 'post');
        
        fov_mask = u_load([vessel_fov_dir vessel_sizes{i_sz} v_files(i_v).name(1:end-12) '_f_mask.mat']);
        fov_mask = padarray(padarray(fov_mask, [pad_t pad_l], 0, 'pre'), [pad_b pad_r], 0, 'post');
        
        contour_struc.outer_edge(:,1) = contour_struc.outer_edge(:,1)+pad_l;
        contour_struc.outer_edge(:,2) = contour_struc.outer_edge(:,2)+pad_t;
        
        contour_struc.inner_edge(:,1) = contour_struc.inner_edge(:,1)+pad_l;
        contour_struc.inner_edge(:,2) = contour_struc.inner_edge(:,2)+pad_t;
        
        contour_struc.vessel_centre(:,1) = contour_struc.vessel_centre(:,1)+pad_l;
        contour_struc.vessel_centre(:,2) = contour_struc.vessel_centre(:,2)+pad_t;
        
        contour_struc.rows = er1 - sr1 + 1;
        contour_struc.cols = ec1 - sc1 + 1;
        
        %Save everything
        display(['Saving copies of image ' num2str(im_number)]);
        im_number = im_number+1;
        save([vessel_mask_dirc vessel_sizes{i_sz} v_files(i_v).name(1:end-12) '_v_mask.mat'], 'vessel_mask');
        save([vessel_cmask_dirc vessel_sizes{i_sz} v_files(i_v).name(1:end-12) '_v_cmask.mat'], 'vessel_centre_mask');
        save([vessel_width_dirc vessel_sizes{i_sz} v_files(i_v).name(1:end-12) '_width.mat'], 'width_map');
        save([vessel_ori_dirc vessel_sizes{i_sz} v_files(i_v).name(1:end-12) '_ori.mat'], 'ori_map');
        save([vessel_fov_dirc vessel_sizes{i_sz} v_files(i_v).name(1:end-12) '_f_mask.mat'], 'fov_mask');      
        save([vessel_image_dirc vessel_sizes{i_sz} v_files(i_v).name(1:end-12) '.mat'], 'vessel_patch');
        
        save([vessel_cont_dirc vessel_sizes{i_sz} v_files(i_v).name(1:end-12) '.mat'], 'contour_struc');
        
        if false && i_v < 10
            figure;
            subplot(1,2,1); imgray(vessel_patch);
            plot(contour_struc.inner_edge(:,1), contour_struc.inner_edge(:,2));
            plot(contour_struc.outer_edge(:,1), contour_struc.outer_edge(:,2));
            plot(contour_struc.vessel_centre(:,1), contour_struc.vessel_centre(:,2));
            plot(contour_struc.vessel_centre(contour_struc.apex_idx,1), contour_struc.vessel_centre(contour_struc.apex_idx,2), 'rx');
            
            subplot(1,2,2); imgray(vessel_mask);
            plot(contour_struc.inner_edge(:,1), contour_struc.inner_edge(:,2));
            plot(contour_struc.outer_edge(:,1), contour_struc.outer_edge(:,2));
            plot(contour_struc.vessel_centre(:,1), contour_struc.vessel_centre(:,2));
            plot(contour_struc.vessel_centre(contour_struc.apex_idx,1), contour_struc.vessel_centre(contour_struc.apex_idx,2), 'rx');
        end
    end
end
%%
%%
dir_list = {'fov_masks', 'images', 'vessel_centre_masks', 'vessel_masks', 'width_maps', 'orientations'};

orig_dir = 'C:\isbe\nailfold\data\rsa_study\cxx\';
new_dir = 'C:\isbe\nailfold\data\rsa_study\cxx2\';
create_folder(new_dir);
for i_dir = 1%2:length(dir_list)
    
    create_folder([new_dir dir_list{i_dir}]);
    
    file_list = dir([orig_dir dir_list{i_dir} '\*mat']);
    
    for i_file = 1:length(file_list)
  
        vessel_im = u_load([orig_dir dir_list{i_dir} '\' file_list(i_file).name]);
        vessel_im = imresize(vessel_im, 0.5, 'lanczos2');
        
        if strcmpi(dir_list{i_dir}(1), 'o')
            new_name = [new_dir dir_list{i_dir} '\' file_list(i_file).name(1:end-4) '.txt'];
            write_complex_txt(vessel_im, new_name);
            
        elseif strcmpi(dir_list{i_dir}(1), 'w')
            new_name = [new_dir dir_list{i_dir} '\' file_list(i_file).name(1:end-4) '.txt'];
            write_array_txt(vessel_im, new_name);
            
        elseif strcmpi(dir_list{i_dir}(1), 'i')
            new_name = [new_dir dir_list{i_dir} '\' file_list(i_file).name(1:end-4) '.png'];
            imwrite(uint8(vessel_im), new_name);
            
        else
            new_name = [new_dir dir_list{i_dir} '\' file_list(i_file).name(1:end-4) '.png'];
            vessel_im = uint8(255*vessel_im);
            imwrite(vessel_im, new_name);
            
        end
        
        %Now take a flipped copy
        vessel_im = fliplr( vessel_im );
        
        if dir_list{i_dir}(1) == 'o'
            vessel_im = complex(real(vessel_im), -imag(vessel_im));
        end      
        
        if strcmpi(dir_list{i_dir}(1), 'o')
            new_name = [new_dir dir_list{i_dir} '\r_' file_list(i_file).name(1:end-4) '.txt'];
            write_complex_txt(vessel_im, new_name);
            
        elseif strcmpi(dir_list{i_dir}(1), 'w')
            new_name = [new_dir dir_list{i_dir} '\r_' file_list(i_file).name(1:end-4) '.txt'];
            write_array_txt(vessel_im, new_name);
            
        elseif strcmpi(dir_list{i_dir}(1), 'i')
            new_name = [new_dir dir_list{i_dir} '\r_' file_list(i_file).name(1:end-4) '.png'];
            imwrite(uint8(vessel_im), new_name);
            
        else
            new_name = [new_dir dir_list{i_dir} '\r_' file_list(i_file).name(1:end-4) '.png'];
            vessel_im = uint8(255*vessel_im);
            imwrite(vessel_im, new_name);
            
        end
        
    end
end
%%
% Write out the training masks and output files in text format to be read in C++
base_dir = 'C:/isbe/nailfold/data/rsa_study/cxx2/';
im_types = {...
        'images/', '.mat', 'images_png/', '.png';...
        'fov_masks/', '_f_mask.mat', 'fov_masks_cxx/', '_f_mask.png';...
        'vessel_masks/', '_v_mask.mat', 'vessel_masks_cxx/', '_v_mask.png';
        'vessel_centre_masks/', '_v_cmask.mat', 'vessel_centre_masks_cxx/', '_vc_mask.png';
        'width_maps/', '_width.mat', 'width_maps_cxx/', '_width.txt';
        'orientations/', '_ori.mat', 'orientations_cxx/', '_ori.txt'};
num_types = size(im_types,1);    
for i_type = 1:num_types
    create_folder([base_dir im_types{i_type,3}]);
end
    
im_list = dir([base_dir 'images\*.mat']);
for i_im = 1:length(im_list);
    im_name = im_list(i_im).name(1:end-4);


    for i_type = 1:num_types
        
        orig_name = [base_dir im_types{i_type,1} im_name im_types{i_type,2}];
        new_name  = [base_dir im_types{i_type,3} im_name im_types{i_type,4}];
        
        vessel_im = u_load(orig_name);
        if strcmpi(im_types{i_type,1}(1), 'o')
            write_complex_txt(vessel_im, new_name);
            
        elseif strcmpi(im_types{i_type,1}(1), 'w')
            save(new_name, 'vessel_im', '-ascii','-double');
            
        elseif strcmpi(im_types{i_type,1}(1), 'i')
            imwrite(uint8(vessel_im), new_name);
            
        else
            vessel_im = uint8(255*vessel_im);
            imwrite(vessel_im, new_name);
            
        end
        
    end
end
%%
base_dir = 'C:/isbe/nailfold/data/rsa_study/cxx2/';
contour_dir = 'C:/isbe/nailfold/data/rsa_study/cxx/vessel_contours/';
new_contour_dir = 'C:/isbe/nailfold/data/rsa_study/cxx2/vessel_contours/';
create_folder(new_contour_dir);

contour_list = dir([contour_dir '*.mat']);

for i_c = 1:length(contour_list);
    
    vessel_contour = u_load([contour_dir contour_list(i_c).name]);
    write_vessel_contour_txt(vessel_contour, [new_contour_dir contour_list(i_c).name(1:end-4) '.txt'], 0.5, 0);
    write_vessel_contour_txt(vessel_contour, [new_contour_dir 'r_' contour_list(i_c).name(1:end-4) '.txt'], 0.5, 1);
end
%%
old_base_dir = 'C:/isbe/nailfold/data/rsa_study/set12g_half/';
new_base_dir = 'C:/isbe/nailfold/data/rsa_study/cxx2/';

transfer_folders = {
    'predictions\detection\rf_classification\296655\';
    'predictions\orientation\rf_regression\296621\';
    'predictions\width\rf_regression\297037\'};

for i_dir = 1:length(transfer_folders)
    
    create_folder([new_base_dir transfer_folders{i_dir}]);
    
    im_list = dir([old_base_dir transfer_folders{i_dir} '*.mat']);
    
    for i_im = 1:length(im_list);
        
        pred_im = u_load([old_base_dir transfer_folders{i_dir} im_list(i_im).name]);
        new_name = [new_base_dir transfer_folders{i_dir} im_list(i_im).name(1:end-4) '.txt'];
        
        if i_dir == 2
            write_complex_txt(pred_im, new_name);
            
        else
            write_array_txt(pred_im, new_name);
            
        end
    end
end
       
%%
old_base_dir = 'C:/isbe/nailfold/data/rsa_study/set12g_half/';
new_base_dir = 'C:/isbe/nailfold/data/rsa_study/cxx2/';

transfer_folders = {
    'predictions\detection\rf_classification\296655\';
    'predictions\orientation\rf_regression\296621\';
    'predictions\width\rf_regression\297037\'};

im_number = 1;
for i_sz = 1:3

    vessel_dir = ['C:\isbe\nailfold\data\rsa_study\apexes\' vessel_sizes{i_sz} '\'];
    contour_dir = ['C:\isbe\nailfold\data\rsa_study\vessel_contours\' vessel_sizes{i_sz} '\'];

    v_files = dir([contour_dir '*contour.mat']);
    for i_v = 1:length(v_files)
        im_name = v_files(i_v).name(1:end-12);
        
        contour_struc = load([contour_dir v_files(i_v).name]);
        vessel_struc = u_load([vessel_dir im_name '.mat']);
        
        full_image = imread(vessel_struc.vessel_properties.nailfold_name);

        [rows cols] = size(full_image);
        [patch_rows patch_cols] = size(vessel_struc.vessel_patch);
        
        sr0 = vessel_struc.vessel_properties.sr;
        er0 = vessel_struc.vessel_properties.sr + patch_rows - 1;
        sc0 = vessel_struc.vessel_properties.sc;
        ec0 = vessel_struc.vessel_properties.sc + patch_cols - 1;
        
        extra_rows = floor((512-patch_rows)/2);
        extra_cols = floor((512-patch_cols)/2);
        if extra_rows <= 0
            sr1 = sr0;
            er1 = er0;
            
        elseif sr0 < extra_rows;
            sr1 = 1;
            er1 = 512;
            
        elseif er0+extra_rows > rows
            er1 = rows;
            sr1 = rows-511;
            
        else
            sr1 = sr0-extra_rows;
            er1 = sr1 + 511;
        end
            
        if extra_cols <= 0
            sc1 = sc0;
            ec1 = ec0;
            
        elseif sc0 < extra_cols;
            sc1 = 1;
            ec1 = 512;
            
        elseif ec0+extra_cols > cols
            ec1 = cols;
            sc1 = cols-511;
            
        else
            sc1 = sc0-extra_cols;
            ec1 = sc1 + 511;
        end
        
        pad_t = floor((sr0 - sr1)/2);
        pad_l = floor((sc0 - sc1)/2);
        
        sized_im_name = [vessel_sizes{i_sz} im_name];
        vessel_mask = imread([new_base_dir 'vessel_masks/' sized_im_name '_v_mask.png']);
        
        vessel_pred = u_load([old_base_dir transfer_folders{1} sized_im_name '_pred.mat']);
        orientation_pred = u_load([old_base_dir transfer_folders{2} sized_im_name '_pred.mat']);
        width_pred = u_load([old_base_dir transfer_folders{3} sized_im_name '_pred.mat']);
        
        pad_b = size(vessel_mask,1) - pad_t - size(vessel_pred,1);
        pad_r = size(vessel_mask,2) - pad_l - size(vessel_pred,2); 
        
        %Now pad everything appropriately
        vessel_pred = padarray(padarray(vessel_pred, [pad_t pad_l], 0, 'pre'), [pad_b pad_r], 0, 'post');
        orientation_pred = padarray(padarray(orientation_pred, [pad_t pad_l], 0, 'pre'), [pad_b pad_r], 0, 'post');
        width_pred = padarray(padarray(width_pred, [pad_t pad_l], 0, 'pre'), [pad_b pad_r], 0, 'post');
        
        write_array_txt(vessel_pred, [new_base_dir transfer_folders{1} sized_im_name '_pred.txt']);
        write_complex_txt(orientation_pred, [new_base_dir transfer_folders{2} sized_im_name '_pred.txt']);
        write_array_txt(width_pred, [new_base_dir transfer_folders{3} sized_im_name '_pred.txt']);
            

        
        if i_v < 10
            figure;
            subplot(1,2,1); imgray(vessel_mask);           
            subplot(1,2,2); imgray(vessel_pred);
            
        end
    end
end
%%
extract_vessel_centres_set(...
    'num_jobs', 1, ...
    'task_id', 1,...
    'data_dir',             [nailfoldroot 'data/rsa_study/cxx2/'],...
    'data_ext',             '*.txt',...
    'prob_dir',             'rf_classification/296655',...
    'ori_dir',              'rf_regression/296621',...
    'width_dir',            'rf_regression/297037',...
    'fov_mask_dir',         '',...
    'centre_dir',           'vessel_centres',...
    'overwrite',            0);
%%
%Write out the forests for preidcting vessel properties and apex location in text format
% So they can be read into the cxx program
complex_rep = 2;
num_angles = 6; 
num_levels = 5;
win_size = 9; %3x3

model_names = {
    'C:\isbe\nailfold\models\vessel\width\rf_regression\297037\predictor.mat'
	'C:\isbe\nailfold\models\vessel\orientation\rf_regression\296621\predictor.mat'
	'C:\isbe\nailfold\models\vessel\detection\rf_classification\296655\predictor.mat'
    'C:\isbe\nailfold\models\vessel\width\rf_regression\675754\predictor.mat'
    'C:\isbe\nailfold\models\vessel\orientation\rf_regression\675752\predictor.mat'
    'C:\isbe\nailfold\models\vessel\detection\rf_classification\675753\predictor.mat'
	'C:\isbe\nailfold\models\apex\classification\frog\rf.mat'
    'C:\isbe\nailfold\models\apex\classification\set12g_half_296655\rf.mat'
    'C:\isbe\nailfold\models\apex\offset_x\set12g_half_296655\rf.mat'
    'C:\isbe\nailfold\models\apex\offset_y\set12g_half_296655\rf.mat'
    'C:\isbe\nailfold\models\apex\rescoring\miccai_all\rf.mat'
    'C:\isbe\nailfold\models\apex\rescoring\corrected_miccai_all\rf.mat'};

tree_roots = {
    'C:\isbe\nailfold\models\vessel\width\rf_regression\'
	'C:\isbe\nailfold\models\vessel\orientation\rf_regression\'
	'C:\isbe\nailfold\models\vessel\detection\rf_classification\'
    'C:\isbe\nailfold\models\vessel\width\rf_regression\'
	'C:\isbe\nailfold\models\vessel\orientation\rf_regression\'
	'C:\isbe\nailfold\models\vessel\detection\rf_classification\'
	''
    ''
    ''
    ''
    ''
    ''};


for i_m = 11%1:length(model_names)
    
    if (ismember(i_m, [1 2 3]))
        translated_dims = translate_feature_dimensions(complex_rep,num_angles,num_levels,win_size); 
    elseif (ismember(i_m, [4 5 6]))
        translated_dims = translate_feature_dimensions(complex_rep,num_angles,num_levels,1);
    else
        translated_dims = -1;
    end

    rf = u_load(model_names{i_m});
    if ~isempty(tree_roots{i_m})
        rf.tree_root = tree_roots{i_m};
    end
    write_forest_txt(rf, [rf.tree_root rf.tree_dir 'txt'], translated_dims);
    create_folder([rf.tree_root rf.tree_dir 'bfs']);
end
%%
base_dir = 'C:/isbe/nailfold/data/rsa_study/master_set/';
im_types = {...
        'images/', '.mat', 'images_png/', '.png';...
        'fov_masks/', '_f_mask.mat', 'fov_masks_cxx/', '_f_mask.png';...
        'predictions\detection\rf_classification\296655\', '_pred.mat', 'predictions\detection\rf_classification\296655_cxx\', '_pred.txt';
        'predictions\orientation\rf_regression\296621\', '_pred.mat', 'predictions\orientation\rf_regression\296621_cxx\', '_pred.txt';
        'predictions\width\rf_regression\297037\', '_pred.mat', 'predictions\width\rf_regression\297037_cxx\', '_pred.txt'};
num_types = size(im_types,1);    
for i_type = 1:num_types
    create_folder([base_dir im_types{i_type,3}]);
end
    
%im_list = dir([base_dir 'images\*.mat']);
load('C:\isbe\nailfold\data\rsa_study\data_lists\miccai_lists.mat', 'miccai_selection');
load('C:\isbe\nailfold\data\rsa_study\data_lists\image_id_data.mat');
im_names = sort(image_id_data.im_names(miccai_selection.test));

for i_im = 1:length(im_names);
    im_name = im_names{i_im}; %im_list(i_im).name(1:end-4);
    display(['converting image ' num2str(i_im)]);

    for i_type = 1:2%num_types
        
        orig_name = [base_dir im_types{i_type,1} im_name im_types{i_type,2}];
        new_name  = [base_dir im_types{i_type,3} im_name im_types{i_type,4}];
        
        if exist(new_name, 'file')
            continue;
        end
        
        vessel_im = u_load(orig_name);
        switch i_type
            case 1
                imwrite(uint8(vessel_im), new_name);
            case 2
                vessel_im = uint8(255*vessel_im);
                imwrite(vessel_im, new_name);
            case {3,5}
                write_array_txt(vessel_im, new_name);
            case 4
                write_complex_txt(vessel_im, new_name);
        end
        
    end
end
%%
