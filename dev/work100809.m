rf_list = dir('M:\asymmetry_project\data\contralateral_rfs\2004_screening\contralateral_images\abnormals\random_forest*.mat');
for ii = 1:length(rf_list)   

    load(['M:\asymmetry_project\data\contralateral_rfs\2004_screening\contralateral_images\abnormals\' rf_list(ii).name]);
    prob1 = random_forest.image1_votes1 ./ random_forest.image1_total_votes;
    prob1(~random_forest.image1_total_votes) = 0;
    prob2 = random_forest.image2_votes1 ./ random_forest.image2_total_votes;
    prob2(~random_forest.image2_total_votes) = 0;
    c1 = min([min(prob1(random_forest.image1_total_votes>0)) min(prob2(random_forest.image2_total_votes>0))]);
    c2 = max([max(prob1(random_forest.image1_total_votes>0)) max(prob2(random_forest.image2_total_votes>0))]);
    figure; 
    subplot(1,2,1); imagesc(prob1); axis image; colormap(jet(256)); caxis([c1 c2]); colorbar('west');
    subplot(1,2,2); imagesc(fliplr(prob2)); axis image; colormap(jet(256)); caxis([c1 c2]);  colorbar('east');
end
%%
%Load a 128x128 smooth patch and add a pretty obvious line
bg = u_load('C:\isbe\asymmetry_project\data\synthetic_backgrounds\smooth128\bg001.mat');
[cls label] = create_gauss_bar(3, 8, 23, 128, 128, 64, 64);
cls = bg+cls;
figure; imagesc(cls); axis image;
%
rf = u_load('C:\isbe\dev\classification\rf\rf_linop_3_5_8\random_forest.mat');
rf.nclasses = rf.n_classes; rf = rmfield(rf, 'n_classes');
prob1 = classify_image_linop('image_in', cls, 'forest', rf, 'use_probs', 0);
prob2 = classify_image_linop('image_in', cls, 'forest', rf, 'use_probs', 1);
%%
mb_combine_rfs(...
    'rf_dir', 'Z:\asymmetry_project\data\line_detection_rfs\191562\',...
    'tree_dir', '191562\',...
    'tree_root', [asymmetryroot 'data/line_detection_rfs/'],...
    'replace_tree_root', 'Z:\asymmetry_project\data\line_detection_rfs\',...
    'copy_trees', 1,...% the optional arguments
    'delete_trees', 0,...
    'save_path', [asymmetryroot 'data\line_detection_rfs\191562\']);
%%
rf = u_load('C:\isbe\asymmetry_project\data\line_detection_rfs\191562\random_forest.mat');
bg = u_load('C:\isbe\asymmetry_project\data\synthetic_backgrounds\smooth128\bg001.mat');
%
[cls label] = create_gauss_bar(3, 16, 23, 128, 128, 64, 64);
cls = bg+cls;
%%
prob1 = classify_image('image_in', cls, 'forest', rf, 'use_probs', 0);
prob2 = classify_image('image_in', cls, 'forest', rf, 'use_probs', 1);

figure; 
subplot(2,2,1:2); imagesc(cls); axis image;
subplot(2,2,3); imagesc(prob1); axis image;
subplot(2,2,4); imagesc(prob2); axis image;
%%
prob1 = 1 - random_forest.image1_votes1 ./ random_forest.image1_total_votes;
prob1(~random_forest.image1_total_votes) = 0;
prob2 = 1 - random_forest.image2_votes1 ./ random_forest.image2_total_votes;
prob2(~random_forest.image2_total_votes) = 0;

figure;
subplot(1,2,1); imagesc(prob1); axis image; colormap(jet(256));
subplot(1,2,2); imagesc(prob2); axis image; colormap(jet(256));
%%
mb_combine_rfs(...
    'rf_dir', 'Z:\asymmetry_project\data\line_detection_rfs\191613\',...
    'tree_dir', '191613\',...
    'tree_root', [asymmetryroot 'data/line_detection_rfs/'],...
    'replace_tree_root', 'Z:\asymmetry_project\data\line_detection_rfs\',...
    'copy_trees', 1,...% the optional arguments
    'delete_trees', 0,...
    'save_path', [asymmetryroot 'data\line_detection_rfs\191613\']);