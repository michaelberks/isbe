load('C:\isbe\nailfold\data\rsa_study\data_lists\image_id_data.mat');
im_count = 0;

im_dir = 'C:\isbe\nailfold\data\rsa_study\images\';
new_dir = 'C:\isbe\cxx\bin32\vxl\external\isbe_apps\nailfold\qtools\ncm_qseries\images\';
categories = cell(10,1);
people_ids = zeros(10,1);
top_older = false(10,1);

for pp = 1:max(image_id_data.people_id)
    
    im_idx1 = ...
        image_id_data.people_id == pp &...
        image_id_data.visit == 1 &...
        image_id_data.digit == 4 &...
        strcmpi(image_id_data.hand, 'L');
    
    im_idx2 = ...
        image_id_data.people_id == pp &...
        image_id_data.visit == 3 &...
        image_id_data.digit == 4 &...
        strcmpi(image_id_data.hand, 'L');
    
    if sum(im_idx1 == 1) && sum(im_idx2 == 1)
        name1 = image_id_data.im_names(im_idx1);
        name2 = image_id_data.im_names(im_idx2);
        
        name1 = [name1{1} '.png'];
        name2 = [name2{1} '.png'];
        
        if exist([im_dir name1], 'file') && exist([im_dir name2], 'file')
            im_count = im_count + 1;
            top_older(im_count) = rand > 0.5;
            
            if top_older(im_count)
                display([name1 ' ' name2]);
            else
                display([name2 ' ' name1]);
            end
            people_ids(im_count) = pp;
            
            cc = image_id_data.category(im_idx1);
            categories(im_count) = cc;
            
            %copyfile([im_dir name1], [new_dir name1]);
            %copyfile([im_dir name2], [new_dir name2]);
            
        end
    end
    
    if im_count == 10
        break;
    end
end
%%
study_dir = 'C:\isbe\nailfold\data\2_year_study\';
[~, ~, raw] = xlsread([study_dir '\data_lists\2yr_data.xls']);
%%
subject_id = raw(2:end,1);
include =  logical(cell2mat(raw(2:end,20)));

%subject_id = subject_id(include);

num_subjects = length(subject_id);
im_dir = [study_dir '\anonymous_png\'];

include_found = false(num_subjects, 1);
for pp = 1:num_subjects
    
    name1 = [subject_id{pp} 'V1LD4X3LrgMosaic.png'];
    name2 = [subject_id{pp} 'V6LD4X3LrgMosaic.png']; 
    
    include_found(pp) =...
        exist([im_dir name1], 'file') && exist([im_dir name2], 'file');    
    
end
%%
subject_id = subject_id(include_found);
num_subjects = length(subject_id);

subject_order = randperm(num_subjects);
subject_id = subject_id(subject_order);
top_older = false(num_subjects,1);

for pp = 1:num_subjects
    
    name1 = [subject_id{pp} 'V1LD4X3LrgMosaic.png'];
    name2 = [subject_id{pp} 'V6LD4X3LrgMosaic.png'];
    
    top_older(pp) = rand > 0.5;
            
    if top_older(pp)
        display([name1 ' ' name2]);
    else
        display([name2 ' ' name1]);
    end
end
save([study_dir '\data_lists\subject_order.mat'], 'subject_order', 'subject_id', 'top_older');
%%
%Loading all the data into an excel spreadsheet
grade_dir = 'N:\musculoskeletal\NCM only studies\2 yr data incl marina mark up\grades2';
xls_filename = [grade_dir '\v1_v6_comparison_a.xls'];
load([study_dir '\data_lists\subject_order.mat'], 'subject_order', 'subject_id', 'top_older');
load_qseries_grades_to_xls(grade_dir, xls_filename, subject_id, top_older);
%%
%--------------------------------------------------------------------------
missing = include & ~include_found;
missing_ids = subject_id(missing);
num_missing = length(missing_ids);
missing_found = false(num_missing,1);
missing_names = cell(num_missing,1);
for pp = 1:num_missing
    
    ns1 = dir([im_dir missing_ids{pp} 'V1*.png']);
    ns2 = dir([im_dir missing_ids{pp} 'V6*.png']); 

    if length(ns1)==1 && length(ns2)==1 && ...
            strcmpi(ns1.name(end-17:end), ns2.name(end-17:end))
        missing_found(pp) = 1;
        missing_names{pp,1} = ns1(1).name;
        missing_names{pp,2} = ns2(1).name;
    end
end
%%
subject_id2 = missing_ids(missing_found);
subject_names = missing_names(missing_found,:);

num_subjects = length(subject_id2);

subject_order2 = randperm(num_subjects);
subject_id2 = subject_id2(subject_order2);
top_older2 = false(num_subjects,1);

for pp = 1:num_subjects
    
    top_older2(pp) = rand > 0.5;
            
    if top_older2(pp)
        display([subject_names{pp,1} ' ' subject_names{pp,2}]);
    else
        display([subject_names{pp,2} ' ' subject_names{pp,1}]);
    end
end
save([study_dir '\data_lists\subject_order2.mat'], 'subject_order2', 'subject_id2', 'top_older2');
    
    
