pair_list = dir('G:\classification_maps\2004_screening\normal_roi_all\image_pair*');

for ii = 1:length(pair_list)
    movefile(...
        ['G:\classification_maps\2004_screening\normal_roi_all\' pair_list(ii).name '\votes_image.mat'],...
        ['G:\classification_maps\2004_screening\normal_roi_all\' pair_list(ii).name '\pair_dataW3L5M0.mat']);
end
%%
pair_list = dir('M:\asymmetry_project\data\contralateral_rfs\2004_screening\contralateral_abnormal_max23July_rfs\image_pair*');

for ii = 1:length(pair_list)
    copyfile(...
        ['M:\asymmetry_project\data\contralateral_rfs\2004_screening\contralateral_abnormal_max23July_rfs\' pair_list(ii).name '\votes_image.mat'],...
        ['G:\classification_maps\2004_screening\abnormal_roi\' pair_list(ii).name '\pair_dataW3L5M1.mat']);
end
%%