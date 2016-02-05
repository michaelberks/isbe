base_dir = 'A:\D_2.2_Datasets\125_375_Single_MLO_Dataset_VOLPARA\';
control_list = dir([base_dir 'CONTROLS\ASSURE_CONTROLS*']);
cancer_list = dir([base_dir 'CANCERS\ASSURE_CANCERS*']);

num_controls = length(control_list);
num_cancers = length(cancer_list);
%
control_im_names = cell(num_controls,1);
control_mask_names = cell(num_controls,1);
for i_case = 1:num_controls
    im_name = dir([base_dir 'CONTROLS\' control_list(i_case).name '\*ML*hint_densityMap*.pgm']);
    mask_name = dir([base_dir 'CONTROLS\' control_list(i_case).name '\*ML*hint_Segmentation.pgm']);
    
    control_im_names{i_case} = [base_dir 'CONTROLS\' control_list(i_case).name '\' im_name(end).name];
    control_mask_names{i_case} = [base_dir 'CONTROLS\' control_list(i_case).name '\' mask_name(end).name];
end

cancer_im_names = cell(num_cancers,1);
cancer_mask_names = cell(num_cancers,1);
for i_case = 1:num_cancers
    im_name = dir([base_dir 'CANCERS\' cancer_list(i_case).name '\*ML*hint_densityMap*.pgm']);
    mask_name = dir([base_dir 'CANCERS\' cancer_list(i_case).name '\*ML*hint_Segmentation.pgm']);
    
    cancer_im_names{i_case} = [base_dir 'CANCERS\' cancer_list(i_case).name '\' im_name(end).name];
    cancer_mask_names{i_case} = [base_dir 'CANCERS\' cancer_list(i_case).name '\' mask_name(end).name];
end
    
case_im_names = [control_im_names; cancer_im_names];
case_mask_names = [control_mask_names; cancer_mask_names];

num_cases = length(case_im_names);
rand_i = randperm(num_cases);

num_train = num_cases*0.8;
num_test = num_cases - num_train;
train_case_idx = rand_i(1:num_train);
test_case_idx = rand_i(num_train+1:end);
%%
fold_idx = load('A:\D_2.2_Datasets\5FCV_125_AgeMatched_indices.mat');

cancer_fold_names = cell(num_cancers,1);%125
control_fold_names = cell(num_controls,1);%375

for i_fold = 1:5
    cancer_idx = fold_idx.cancer_fold_indx(i_fold,:);
    
    for i_can = 1:25
        
        %Assign new name for this cancer
        cancer_fold_names{cancer_idx(i_can)} = ['fold_' num2str(i_fold) '_cancer_' zerostr(i_can, 3)];
        
        %Find matching controls
        control_idx = find(fold_idx.Matched_Controls_Indices==cancer_idx(i_can));
        
        if length(control_idx) ~= 3
            display(['Wrong number (' num2str(length(control_idx)) ') matched controls found for cancer #' num2str(cancer_idx(i_can))]);
        else
            for i_con = 1:3
                control_fold_names{control_idx(i_con)} = ['fold_' num2str(i_fold) '_control_' zerostr(i_can, 3) '_' zerostr(i_con, 3)];
            end
        end
    end
end
%%
case_fold_names = [control_fold_names; cancer_fold_names];

image_dir = 'c:\isbe\density\assure\375_125\images\';
full_mask_dir = 'c:\isbe\density\assure\375_125\full_masks\';
compressed_mask_dir = 'c:\isbe\density\assure\375_125\compressed_masks\';
label_mask_dir = 'c:\isbe\density\assure\375_125\label_masks\';

create_folder(image_dir);
create_folder(full_mask_dir);
create_folder(compressed_mask_dir);
create_folder(label_mask_dir);

for i_case = 1:num_cases
    
    %Copy mammogram and mask
    mam_vdm = imread(case_im_names{i_case});
    full_mask = imread(case_mask_names{i_case});
    compressed_mask = full_mask == 4;
    if i_case <= num_controls
        label_mask = false(size(mam_vdm));
    else
        label_mask = true(size(mam_vdm));
    end
    
    save([image_dir case_fold_names{i_case} '.mat'], 'mam_vdm');
    save([full_mask_dir case_fold_names{i_case} '_f_mask.mat'], 'full_mask');
    save([compressed_mask_dir case_fold_names{i_case} '_c_mask.mat'], 'compressed_mask');
    save([label_mask_dir case_fold_names{i_case} '_l_mask.mat'], 'label_mask');
    
end
%%
data_list_dir = 'C:\isbe\density\assure\375_125\data_lists\';
create_folder(data_list_dir);
save([data_list_dir 'case_fold_names.mat'], 'case_fold_names', 'control_fold_names', 'case_fold_names');