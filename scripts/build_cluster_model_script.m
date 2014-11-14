%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare some data for clustering

data_type = 'mc';
m_list = dir('C:/isbe/dev/background/images/', data_type, '/', data_type, '*');

num_images = length(m_list);

%Get mammogram regions from saved annotations and save as .bmp in images
%directory
for ii = 1:num_images
    load(['C:/isbe/dev/background/images/', data_type, '/', m_list(ii).name]);
    figure; imagesc(mass.mass_ROI); colormap(gray(256)); axis image;
    imwrite(mass.mass_ROI, ['C:/isbe/dev/background/images/', data_type, '/', data_type, '',...
        zerostr(ii,3), '.bmp']);
    clear mass;
end
%%
% Half the size of the patches and save in small images directory
% Also useful to calculate the total points
total_pts = 0;
for ii = 1:num_images
    im1 = imread(['C:/isbe/dev/background/images/', data_type, '/', data_type, '',...
        zerostr(ii,3), '.bmp']);
    p_small = imresize(im1, 0.5, 'bicubic');
    imwrite(p_small, ['C:/isbe/dev/background/images_small/', data_type, '/', data_type, '',...
        zerostr(ii,3), '.bmp']);
    
    total_pts = numel(p_small) + total_pts;
    clear im1 p_small;
end

%%
% Show the small patches
for ii = 1:10
    im1 = imread(['C:/isbe/dev/background/images_small/', data_type, '/', data_type, '',...
        zerostr(ii,3), '.bmp']);
    figure; imagesc(im1); colormap(gray(256)); axis image;
    clear mass;
end
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build the args structures to build cluster models

data_type = 'pectoral_8';

% Arguments to build indexes for the divided sets
cluster_idx_args.NumberSets = 1;
cluster_idx_args.IdxDir = [mberksroot, 'background/idx/', data_type, '/'];
cluster_idx_args.NextIdxFunctionArgs.ImageDir = [mberksroot, 'background/images/', data_type, '/'];
cluster_idx_args.NextIdxFunctionArgs.TempStorageDir = [mberksroot, 'background/temp/'];
cluster_idx_args.NextIdxFunctionArgs.WindowSize = 15;
cluster_idx_args.NextIdxFunctionArgs.MaxMemory = 256;
cluster_idx_args.NextIdxFunction = 'mb_next_idx_for_clustering';
%%
% Now do the indexing
mb_save_cluster_idx(cluster_idx_args);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Final and divide cluster args a built in the cluster scripts, but can be
% useful to have them here for testing/debugging

final_cluster_args.ClusteringFunctionArgs.k = 40;
final_cluster_args.ClusteringFunctionArgs.MaxIter = 2;
final_cluster_args.FinalFile = [mberksroot, 'background/results/', data_type, '_model'];
final_cluster_args.ResultsDir = [mberksroot, 'background/results/', data_type, '/'];
final_cluster_args.MaxFinalMemory = 256;

%%
divide_cluster_args.IdxFile = [mberksroot, 'background/idx/', data_type, '/cluster_idx001'];
divide_cluster_args.ResultFile = [mberksroot, 'background/results/', data_type, 'clustering_result'];
divide_cluster_args.ClusteringFunctionArgs.k = 40;
divide_cluster_args.ClusteringFunctionArgs.MaxIter = 2;
divide_cluster_args.ClusteringFunctionArgs.PropRepresentativePoints = 0.2;
divide_cluster_args.ImageDir = [mberksroot, '/background/images/', data_type, '/'];
divide_cluster_args.WindowSize = 15;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now we've got some models, let's do some synthesising!
%%
N = 1 %num of background so synthesize
image_dir = 'C:/isbe/dev/background/images_small/normal/'; %Which patches to use
model_path = 'C:/isbe/dev/background/results/mc_model'; %which model to use
syn_dir = 'C:/isbe/dev/background/syn/mc/'; %where to save the synthesised images

m_list = dir([image_dir, '*.bmp']);

for ii = 1:N
    bg_patch = imread([image_dir, m_list(ii).name]);
    [r_sub c_sub] = ind2sub(size(bg_patch), 1:numel(bg_patch));
    
    circ_idx = (c_sub - 150).^2 + (r_sub - 150).^2 < 50.^2;
    filled_image = ones(size(bg_patch));
    filled_image(circ_idx) = 0;
    
    
    [synthesised_image] = mb_gmm_tex_synthesis(...
        'PathToTextureGMM', model_path,...
        'SeededImage', bg_patch,...
        'FilledImage', filled_image,...
        'SaveFrequency', 1000);
    
    save([syn_dir, 'syn_', zerostr(ii, 3)], 'synthesised_image');
    figure; 
    subplot(1,2,1); imagesc(bg_patch); colormap(gray(256)); axis image;
    subplot(1,2,2); imagesc(synthesised_image); colormap(gray(256)); axis image;
    
    clear bg_patch r_sub c_sub circ_idx filled_image;
end
%%
data_type = 'normal_2';
final_cluster_args.ClusteringFunctionArgs.k = 40;
final_cluster_args.ClusteringFunctionArgs.MaxIter = 2;
final_cluster_args.FinalFile = [mberksroot, 'background/results/', data_type, '_model'];
final_cluster_args.MaxFinalMemory = 128;
final_cluster_args.PyramidsDir = [mberksroot, 'background/pyramid/', data_type, '/'];
final_cluster_args.Level = 2;
final_cluster_args.Orientation = 1;