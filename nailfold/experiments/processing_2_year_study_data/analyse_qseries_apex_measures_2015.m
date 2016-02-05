
% At some point double check the creation of the spreadsheet for FI's
% grades. For now, take it as given from the spreadsheet below
%%
%Create the spreadsheet for AH's grades
study_dir = 'C:\isbe\nailfold\data\2_year_study\';
observer = {'aherrick', 'fignoli'};
for i_ob = 2
    grade_dir = [study_dir 'grades\' observer{i_ob} '\'];
    xls_filename = [grade_dir '\' observer{i_ob} '_v1_v6_comparison.xls'];
    load([study_dir '\data_lists\' observer{i_ob} '_subject_order.mat'], 'subject_order', 'subject_id', 'top_older');
    load_qseries_grades_to_xls(grade_dir, xls_filename, subject_id, top_older);
end
%%
num_images = 145;
im_class = zeros(num_images,1);
top_older = false(num_images, 2);
is_gradeable = false(num_images, 2);
is_normal = false(num_images, 2);
progressed = zeros(num_images, 2);
reasons = false(num_images, 13, 2);
    
for i_ob = 1:2
    grade_dir = [study_dir 'grades\' observer{i_ob} '\'];
    [~,~,grades_data] = xlsread([grade_dir '\' observer{i_ob} '_v1_v6_comparison.xls']);
 
    if i_ob == 1
        subject_ids = grades_data(2:end,1);
        im_class(strncmpi(subject_ids, 'c', 1)) = 1;
        im_class(strncmpi(subject_ids, 'l', 1)) = 2;
        im_class(strncmpi(subject_ids, 'd', 1)) = 3;
        im_class(strncmpi(subject_ids, 'p', 1)) = 4;
        im_class(strncmpi(subject_ids, 'u', 1)) = 5;
    else
        if any(~strcmpi(subject_ids, grades_data(2:end,1)))
            display('Subject IDs do not match. Houston, we have a problem.');
        end
    end
       
    top_older(:,i_ob) = strcmpi(grades_data(2:end,2), 'V1');
    is_gradeable(:,i_ob) = cell2mat(grades_data(2:end,3))>0;
    is_normal(:,i_ob) = cell2mat(grades_data(2:end,4))==0;
    progressed(:,i_ob) = cell2mat(grades_data(2:end,5));
    reasons(:,:,i_ob) = cell2mat(grades_data(2:end,6:18))>0;    
end
%%
sum(is_normal(is_gradeable(:,1),1) & is_normal(is_gradeable(:,1),2))
sum(~is_normal(is_gradeable(:,1),1) & is_normal(is_gradeable(:,1),2))
sum(is_normal(is_gradeable(:,1),1) & ~is_normal(is_gradeable(:,1),2))
sum(~is_normal(is_gradeable(:,1),1) & ~is_normal(is_gradeable(:,1),2))
im_classes = {'c', 'l', 'd', 'p', 'u'};

for i_c = 1:5
    display(['Image class = ' im_classes{i_c}]);
    sum(is_normal(is_gradeable(:,1),1) & is_normal(is_gradeable(:,1),2) & im_class(is_gradeable(:,1))==i_c)
    sum(~is_normal(is_gradeable(:,1),1) & is_normal(is_gradeable(:,1),2) & im_class(is_gradeable(:,1))==i_c)
    sum(is_normal(is_gradeable(:,1),1) & ~is_normal(is_gradeable(:,1),2) & im_class(is_gradeable(:,1))==i_c)
    sum(~is_normal(is_gradeable(:,1),1) & ~is_normal(is_gradeable(:,1),2) & im_class(is_gradeable(:,1))==i_c)
end
%%
idx  = is_gradeable(:,1) & ~is_normal(:,1) & ~is_normal(:,2);
prog_cooc = full(sparse(progressed(idx,2)+3,progressed(idx,1)+3,1,5,5))
%%
progression_reasons = {
    'A general increase in capillary size'
    'An increase in the number of giant capillaries'
    'A decrease in the number of giant capillaries'
    'An increase in the number of haemhorrages'
    'An increase in the number of angiogenic capillaries'
    'An increase in the derangement of capillaries'
    'A decrease in the number of capillaries or an increase in avascular areas'
    'Increase in halo effect'
    'Decrease in halo effect'
    'Other changes related to capillary size'
    'Other changes related to capillary shape'
    'An decrease in the number of angiogenic capillaries'
    'An decrease in the derangement of capillaries'};
%%
idx = find((progressed(:,1) + progressed(:,2))>=3 & is_gradeable(:,1));
for i_im = idx(:)'
    display(subject_ids(i_im));
    for i_ob = 1:2
        display(['Reasons given by O' num2str(i_ob)]);
        progression_reasons(reasons(i_im,:,i_ob))
    end
    display('***************************');
end
%%
idx = find((progressed(:,1) + progressed(:,2))<=-3 & is_gradeable(:,1));
for i_im = idx(:)'
    display(subject_ids(i_im));
    for i_ob = 1:2
        display(['Reasons given by O' num2str(i_ob)]);
        progression_reasons(reasons(i_im,:,i_ob))
    end
    display('***************************');
end
%%
load('C:\isbe\nailfold\data\2_year_study\results\people_stats.mat');
people_idx = zeros(num_images,1);
for i_im = 1:num_images
    people_idx(i_im) = find(strcmpi(people_stats.people_ids_str, subject_ids{i_im}));
end
%%
features = {'mean_weighted_width', 'max_mean_width', 'num_giant_vessels', 'mean_orientation_entropy', 'vessel_density1'};
feature_headers = {'C = 22', 'L = 67', 'D = 12', 'P = 8', 'U = 12', 'All = 121'};
feature_cols = {'Mean'; 'Mean (y2-y0)'; 'Mean |y2-y0|'; 'S.D.'; 'S.D. y0'; 'S.D. y2'};
feature_mat = cell(7,7);
feature_mat(1,2:7) = feature_headers;
feature_mat(2:7,1) = feature_cols;

for i_f = 1:5
    X = squeeze(sum(sum(people_stats.(features{i_f}),3),2));
    X = X(:,[1 6]);
    
    for i_c = 1:6
        if i_c < 6
            idx = im_class==i_c & is_gradeable(:,1);
        else
            idx = is_gradeable(:,1);
        end
        idx = people_idx(idx);
        X_all = X(idx,:);
        feature_mat{2, i_c+1} = nanmean(X_all(:));
        feature_mat{3, i_c+1} = nanmean(X(idx,1)-X(idx,2));
        feature_mat{4, i_c+1} = nanmean(abs(X(idx,1)-X(idx,2)));
        feature_mat{5, i_c+1} = nanstd(X_all(:));
        feature_mat{6, i_c+1} = nanstd(X(idx,1));
        feature_mat{7, i_c+1} = nanstd(X(idx,2));
    end
    display(['Stats for feature ' features{i_f}]);
    display(feature_mat);
    xlswrite('C:\isbe\nailfold\data\2_year_study\results\auto_avg_by_grp.xls', feature_mat, features{i_f});
end
%%
features = {'mean_weighted_width', 'max_mean_width', 'num_giant_vessels', 'mean_orientation_entropy', 'vessel_density1'};
feature_headers = {'-2', '-1', '0', '1', '2'};
feature_mat = cell(6,6);
feature_mat(1,2:6) = feature_headers;
feature_mat(2:6,1) = feature_headers';

for i_f = 1:5
    X = squeeze(sum(sum(people_stats.(features{i_f}),3),2));
    X = X(:,[1 6]);
    
    for i_p1 = -2:2
        for i_p2 = -2:2
            
            i_c = i_p1+4;
            i_r = i_p2+4;
            idx = people_idx(is_gradeable(:,1) & progressed(:,1)==i_p1 & progressed(:,2)==i_p2);
            feature_mat{i_c, i_r} = nanmean(X(idx,2)-X(idx,1));
        end
    end

    display(['Stats for feature ' features{i_f}]);
    display(feature_mat);
    xlswrite('C:\isbe\nailfold\data\2_year_study\results\auto_avg_by_prg.xls', feature_mat, features{i_f});
end    
%%
feature_mat = cell(14,7,3);
feature_headers = {'O1', 'O2', 'O1&O2'};

progressed(isnan(progressed))= 0;
for i_type = 1:3

    feature_mat(1,2,i_type) = feature_headers(i_type);
    feature_mat(1,3:7,i_type) = features;
    feature_mat(2:end,1,i_type) = progression_reasons;

    for i_f = 1:5
        X = squeeze(sum(sum(people_stats.(features{i_f}),3),2));
        X = X(:,[1 6]);

        for i_r = 1:13
            
            switch i_type
                case 1
                    idx = reasons(:,i_r,1) & is_gradeable(:,1);
                    progression_i = sign(progressed(:,1));
                    
                case 2
                    idx = reasons(:,i_r,2) & is_gradeable(:,1);
                    progression_i = sign(progressed(:,2));
                    
                case 3
                    idx = reasons(:,i_r,1) & reasons(:,i_r,2) & is_gradeable(:,1);
                    progression_1 = sign(progressed(:,1));
                    progression_2 = sign(progressed(:,2));
                    idx = idx & (progression_1 == progression_2);
                    progression_i = progression_1;
            end

            feature_mat{i_r+1,2,i_type} = sum(idx);
            
            if feature_mat{i_r+1,2,i_type}
                idx1 = idx & progression_i < 0;
                idx2 = idx & progression_i > 0;
                idx1 = people_idx(idx1);
                idx2 = people_idx(idx2);
                x1 = [X(idx1,1); X(idx2,2)];
                x2 = [X(idx1,2); X(idx2,1)];
                feature_mat{i_r+1,i_f+2,i_type} = mean(x1-x2);
            end
        end
    end
    display(feature_mat(:,:,i_type));
    xlswrite('C:\isbe\nailfold\data\2_year_study\results\auto_avg_by_reason.xls', feature_mat(:,:,i_type), i_type);
end
        
        
        
        
   
       

%%
%Francesca marked 46 images as normal
% 
class_names = {'control', 'lSSc', 'dSSc', 'PRP', 'UCTD'};
display([num2str(sum(is_normal)) ' of ' num2str(num_images) ' images marked as normal']);
for i_class = 1:5
    display([num2str(sum(is_normal & im_class==i_class)) ' of ' num2str(sum(im_class==i_class))...
        ' ' class_names{i_class} ' images marked as normal']);
end
display('***');

num_progressed = sum(abs(progressed(~is_normal))>0);

display([num2str(num_progressed) ' of ' num2str(num_images - sum(is_normal)) ' images marked as showing signs of progression']);
for i_class = 1:5
    display([num2str(sum(abs(progressed(~is_normal & im_class==i_class))>0)) ' of '...
        num2str(sum(~is_normal & im_class==i_class)) ' ' class_names{i_class} ' images marked as progressed']);
end
display('***');
%%
num_fwd_progressed = sum(progressed(~is_normal)>0);

display([num2str(num_fwd_progressed) ' of ' num2str(num_progressed) ' images marked as 24 months progressed']);
for i_class = 1:5
    display([num2str(sum(progressed(~is_normal & im_class==i_class)>0)) ' of ' ...
        num2str(sum(abs(progressed(~is_normal & im_class==i_class))>0)) ' ' class_names{i_class} ' images marked as 24 months progressed']);
end
display('***');

num_rev_progressed = sum(progressed(~is_normal)>0);

display([num2str(num_rev_progressed) ' of ' num2str(num_progressed) ' images marked as baseline progressed']);
for i_class = 1:5
    display([num2str(sum(progressed(~is_normal & im_class==i_class)>0)) ' of ' ...
        num2str(sum(abs(progressed(~is_normal & im_class==i_class))>0)) ' ' class_names{i_class} ' images marked as baseline progressed']);
end
display('***');


