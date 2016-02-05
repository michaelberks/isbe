base_dir = 'A:\PROCAS_CASE_CONTROL_SET\';
control_subdir = 'VolparaOutput952Controls\';
cancer_subdir = 'VolparaOutput318Cancer\';
control_subsubdir = 'PROCAS_CC_CONTROLS_RAW';
cancer_subsubdir = 'PROCAS_CASE_CONTROL';
control_list = dir([base_dir control_subdir control_subsubdir '*']);
cancer_list = dir([base_dir cancer_subdir cancer_subsubdir '*']);

num_controls = length(control_list);
num_cancers = length(cancer_list);
%%
control_im_names = cell(num_controls,1);
control_mask_names = cell(num_controls,1);
for i_case = 1:num_controls
    im_name = dir([base_dir control_subdir control_list(i_case).name '\*ML*hint_densityMap*.pgm']);
    mask_name = dir([base_dir control_subdir control_list(i_case).name '\*ML*hint_Segmentation.pgm']);
    
    if length(im_name) >= 1
        control_im_names{i_case} = [base_dir control_subdir control_list(i_case).name '\' im_name(end).name];
        control_mask_names{i_case} = [base_dir control_subdir control_list(i_case).name '\' mask_name(end).name];
    end
end

cancer_im_names = cell(num_cancers,1);
cancer_mask_names = cell(num_cancers,1);
for i_case = 1:num_cancers
    im_name = dir([base_dir cancer_subdir cancer_list(i_case).name '\*ML*hint_densityMap*.pgm']);
    mask_name = dir([base_dir cancer_subdir cancer_list(i_case).name '\*ML*hint_Segmentation.pgm']);

    if length(im_name) >= 1
        cancer_im_names{i_case} = [base_dir cancer_subdir cancer_list(i_case).name '\' im_name(end).name];
        cancer_mask_names{i_case} = [base_dir cancer_subdir cancer_list(i_case).name '\' mask_name(end).name];
    end
end
    
case_im_names = [control_im_names; cancer_im_names];
case_mask_names = [control_mask_names; cancer_mask_names];

num_cases = length(case_im_names);
%%
image_dir = 'c:\isbe\density\assure\procas318\images\';
full_mask_dir = 'c:\isbe\density\assure\procas318\full_masks\';
compressed_mask_dir = 'c:\isbe\density\assure\procas318\compressed_masks\';
label_mask_dir = 'c:\isbe\density\assure\procas318\label_masks\';
datalist_dir = 'c:\isbe\density\assure\procas318\data_lists\';

create_folder(image_dir);
create_folder(full_mask_dir);
create_folder(compressed_mask_dir);
create_folder(label_mask_dir);
create_folder(datalist_dir);
%%
i_control = 1;
i_cancer = 319;

im_names = cell(num_cases,1);
for i_case = 1:num_cases
    
    display(['Saving image ' num2str(i_case)]);

    if isempty(case_im_names{i_case})
        continue;
    end

    %Copy mammogram and mask
    mam_vdm = imread(case_im_names{i_case});
    full_mask = imread(case_mask_names{i_case});
    compressed_mask = full_mask > 0; %== 4
    if i_case <= num_controls
        label_mask = false(size(mam_vdm));
    else
        label_mask = true(size(mam_vdm));
    end
    
    if i_case <= num_controls
        im_name = ['control' zerostr(i_control, 4)];
        i_control = i_control + 1;
    else
        im_name = ['cancer' zerostr(i_cancer, 4)];
        i_cancer = i_cancer + 1;
    end
    im_names{i_case} = im_name;
    
    save([image_dir im_name '.mat'], 'mam_vdm');
    save([full_mask_dir im_name '_f_mask.mat'], 'full_mask');
    save([compressed_mask_dir im_name '_c_mask.mat'], 'compressed_mask');
    save([label_mask_dir im_name '_l_mask.mat'], 'label_mask');
    
end
%%
for i_case = 1:num_cases
    display(['Saving image ' num2str(i_case)]);
    if isempty(case_im_names{i_case})
        continue;
    end
    mam_vdm = u_load([image_dir im_names{i_case} '.mat']);
    mam_vdm = (double(mam_vdm)./ 65535 - 0.5) .* 2.0;
    
    delete([image_dir im_names{i_case} '.mat']);
    save([image_dir im_names{i_case} '.mat'], 'mam_vdm');
end
%%
save([datalist_dir 'case_names.mat'], 'case_im_names', 'case_mask_names', 'im_names');
%%
score_dir = {...
    'A:\PROCAS_CASE_CONTROL_SET_SCORES_VDM\CONTROLS_VDM_SCORES\'
    'A:\PROCAS_CASE_CONTROL_SET_SCORES_VDM\CANCERS_VDM_SCORES\'};
scores_list{1} = dir([score_dir{1} '*.mat']);
scores_list{2} = dir([score_dir{2} '*.mat']);
%
case_score_names = cell(num_cases,1);
for i_case = 1:num_cases
    
    display(['Matching image ' num2str(i_case)]);
    if isempty(case_im_names{i_case})
        continue;
    end
    if i_case <= num_controls
        use_dir = 1;
    else
        use_dir = 2;
    end
    
    [~,file_name,~] = fileparts(case_im_names{i_case});
    
    for i_score = 1:length(scores_list{use_dir})
        if strncmpi(file_name, scores_list{use_dir}(i_score).name, 45)
            case_score_names{i_case} = [score_dir{use_dir} scores_list{use_dir}(i_score).name];
            scores_list{use_dir}(i_score) = [];
            break;
        end
    end
end
save([datalist_dir 'case_names.mat'], 'case_score_names', '-append');
%%
        