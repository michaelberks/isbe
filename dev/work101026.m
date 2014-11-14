r_list = dir('Z:\asymmetry_project\data\relevance_maps\2004_screening\normals\*.mat');

for ii = 181:188
    p1 = u_load(['Z:\asymmetry_project\data\relevance_maps\2004_screening\normals\' r_list(ii).name]);
    
    figure; imagesc(p1); axis image; colormap(jet(256));
    
    clear p1;
end
pack;
%% These images need inverting
% They must have run with the correct code - probably one node sync'd
% earlier than the others on hydra
% Once inverted, the radial maps can be re-run
invert_idx = [57:60 117 118 133 134 169 170 187 188];
for ii = 1:length(invert_idx)
    load(['Z:\asymmetry_project\data\relevance_maps\2004_screening\normals\' r_list(invert_idx(ii)).name]);
    
    figure;
    subplot(1,2,1); imagesc(probability_image); axis image; colormap(jet(256));
    
    probability_image = 1 - probability_image;
    subplot(1,2,2); imagesc(probability_image); axis image;
    
    save(['Z:\asymmetry_project\data\relevance_maps\2004_screening\normals\' r_list(invert_idx(ii)).name], 'probability_image');
    save(['Z:\asymmetry_project\data\mammograms\2004_screening\normals\results\roi_01\' r_list(invert_idx(ii)).name], 'probability_image');
end
%% These boys are bad...
% We need to re-run: line detection, orientation regression and abnormality
% relevance. We can then re-run the radial maps
bad_idx = [41 42 45 46 49 50 61 62 77 78 82 85 173 174];
mammo_names = get_mammo_info(r_list);

for ii = 1:length(bad_idx)
    mammo = u_load(['M:\asymmetry_project\data\mammograms\2004_screening\normals\' mammo_names{bad_idx(ii)} '.mat']);
    mask = u_load(['C:\isbe\asymmetry_project\data\masks\2004_screening\normals\' mammo_names{bad_idx(ii)} '_mask.mat']);
    
    figure;
    subplot(1,2,1); imagesc(mammo); axis image;
    subplot(1,2,2); imagesc(mask); axis image;
    clear mammo mask;
end
%% These images need rotating through 180 degrees
for ii = 1:12
    load(['C:\isbe\asymmetry_project\data\mammograms\2004_screening\normals\' mammo_names{bad_idx(ii)} '.mat']);
    mammo_small = rot90(mammo_small, 2);
    save(['C:\isbe\asymmetry_project\data\mammograms\2004_screening\normals\' mammo_names{bad_idx(ii)} '.mat'], 'mammo_small');
    save(['Z:\asymmetry_project\data\mammograms\2004_screening\normals\' mammo_names{bad_idx(ii)} '.mat'], 'mammo_small');
    save(['M:\asymmetry_project\data\mammograms\2004_screening\normals\' mammo_names{bad_idx(ii)} '.mat'], 'mammo_small');
    clear mammo_small;
end
%% These images need re-segmenting
mam_list = dir('C:\isbe\mammograms\new_CAD\BMP_2004_Normals\044R*.bmp');
segment_breast_batch('C:\isbe\mammograms\new_CAD\BMP_2004_Normals\',...
    'C:\isbe\mammograms\new_CAD\BMP_2004_Normals\segmentations\',...
    'bmp', mam_list, 1);
%%
mammo_dir = [asymmetryroot 'data/mammograms/2004_screening/normals/']; %
seg_dir = [asymmetryroot 'data/segmentations/2004_screening/normals/']; %
mask_dir = [asymmetryroot 'data/masks/2004_screening/normals/']; %

for ii = 1:2;
    mammo = u_load([mammo_dir mam_list(ii).name(1:end-4) '.mat']);
    seg = u_load([seg_dir mam_list(ii).name(1:end-4) '_segmentation.mat']);

    %Resize segmentations
    seg.breast_border = segment_breast_resize(size(mammo), seg);    

    %create masks of breast region for each mammograms
    mask = roipoly(mammo, seg.breast_border(:,1), seg.breast_border(:,2));

    save([mask_dir mam_list(ii).name(1:end-4) '_mask.mat'], 'mask');
end
%%
%--------------------------------------------------------------------------
%**************************************************************************
%--------------------------------------------------------------------------
%% Now check the abnormals... Phew, they're alright!
r_list = dir('Z:\asymmetry_project\data\relevance_maps\2004_screening\abnormals\*.mat');

for ii = 10*29 + (1:10)
    p1 = u_load(['Z:\asymmetry_project\data\relevance_maps\2004_screening\abnormals\' r_list(ii).name]);
    
    figure; imagesc(p1); axis image; colormap(jet(256));
    title(get_mammo_info(r_list(ii).name));
    clear p1;
end
pack;