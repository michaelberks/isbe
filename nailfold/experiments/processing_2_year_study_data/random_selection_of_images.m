function [image_names, subject_order, subject_id, top_older] = random_selection_of_images(varargin)

args = u_packargs(varargin,... % the user's input
    0, ... % non-strict mode
    'study_dir',	'C:\isbe\nailfold\data\2_year_study\',...
    'user_name',    'mberks',...                
    'plot', 0);


%%
study_dir = args.study_dir;
[~, ~, raw] = xlsread([study_dir 'data_lists\2yr_data.xls']);
%%
subject_id = raw(2:end,1);
include =  logical(cell2mat(raw(2:end,20)));

%subject_id = subject_id(include);

num_subjects = length(subject_id);
im_dir = [study_dir '\anonymous_png\'];

include_found = false(num_subjects, 1);
image_names = cell(num_subjects, 2);
for pp = 1:num_subjects
    
    if ~include(pp)
        continue;
    end
    
    name1 = [subject_id{pp} 'V1LD4X3LrgMosaic.png'];
    name2 = [subject_id{pp} 'V6LD4X3LrgMosaic.png']; 
    
    include_found(pp) =...
        exist([im_dir name1], 'file') && exist([im_dir name2], 'file');
    
    if ~include_found(pp)
        ns1 = dir([im_dir subject_id{pp} 'V1*.png']);
        ns2 = dir([im_dir subject_id{pp} 'V6*.png']); 

        if length(ns1)==1 && length(ns2)==1 && ...
            strcmpi(ns1.name(end-17:end), ns2.name(end-17:end));
        
            include_found(pp) = 1;
            name1 = ns1(1).name;
            name2 = ns2(1).name;
        end
    end
    if include_found(pp)
        image_names{pp, 1} = name1;
        image_names{pp, 2} = name2;
    else
        display(['Images not found for subject ' subject_id{pp}]);
    end
end
%%
if ~exist([study_dir '\data_lists\version.mat'], 'file');
    version = 0;
else
    version = u_load([study_dir '\data_lists\version.mat']);
    version = version + 1;
end

subject_id = subject_id(include_found);
image_names = image_names(include_found,:);
num_subjects = length(subject_id);

subject_order = randperm(num_subjects);
subject_id = subject_id(subject_order); %#ok
image_names = image_names(subject_order,:);
top_older = false(num_subjects,1);

im_list_path = [study_dir '\data_lists\' args.user_name '_images(' num2str(version) ').txt'];   
fid = fopen(im_list_path, 'wt');

for pp = 1:num_subjects
    
    name1 = image_names{pp,1};
    name2 = image_names{pp,2};
    
    top_older(pp) = rand > 0.5;
            
    if top_older(pp)
        display([name1 ' ' name2]);
        fprintf(fid, '%s %s\n', name1, name2);
    else
        display([name2 ' ' name1]);
        fprintf(fid, '%s %s\n', name2, name1);
    end    
end
fclose(fid);

save([study_dir '\data_lists\version.mat'], 'version');
save([study_dir '\data_lists\subject_order( ' num2str(version) ').mat'], 'image_names', 'subject_order', 'subject_id', 'top_older');

    
    
