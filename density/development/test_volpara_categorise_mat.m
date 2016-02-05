load('A:\PROCAS_ALL_VOLPARA_SPREADSHEETS\Volpara_Final_ALL.mat');
%%
num_images = size(Volpara_final_all, 1)-1;

labels = cell(num_images,1);
subject_ids = cell(num_images,1);
for i_im = 1:num_images
    if ~isempty(strfind(Volpara_final_all{i_im+1,9}, 'LCC'))
        labels{i_im} = 'LCC';
    elseif ~isempty(strfind(Volpara_final_all{i_im+1,9}, 'LML'))
        labels{i_im} = 'LML';
    elseif ~isempty(strfind(Volpara_final_all{i_im+1,9}, 'RCC'))
        labels{i_im} = 'RCC';
    elseif ~isempty(strfind(Volpara_final_all{i_im+1,9}, 'RML'))
        labels{i_im} = 'RML';
    else
        display(['No matching view found for image ' num2str(i_im) ': ' Volpara_final_all{i_im+1,9}]);
    end
    subject_ids{i_im} = Volpara_final_all{i_im+1,9}(12:16);
end
densities = cell2mat(Volpara_final_all(2:end,7));
%%
[unique_case_ids, vdg4, vdg5, vbd_mean, vbd_lr_max,...
    vbd_scores_by_case, views_present, outliers] = ...
    volpara_categorise_vbd_batch(densities, subject_ids, labels);



