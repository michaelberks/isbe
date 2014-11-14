data_type = {'abnormals', 'normals'};

for jj = 1:2

    mammo_dir = [asymmetryroot 'data/mammograms/2004_screening/' data_type{jj} '/mat/']; %
    seg_dir = [asymmetryroot 'data/segmentations/2004_screening/' data_type{jj} '/']; %
    mask_dir = [asymmetryroot 'data/masks/2004_screening/' data_type{jj} '/']; %
    mkdir(mask_dir);

    %Get list of mammogrmas
    m_list = dir([mammo_dir '*.mat']);

    for ii = 1:length(m_list);
        mammo = u_load([mammo_dir m_list(ii).name]);
        seg = u_load([seg_dir m_list(ii).name(1:end-4) '_segmentation.mat']);

        %Resize segmentations
        seg.breast_border = segment_breast_resize(size(mammo), seg);    

        %create masks of breast region for each mammograms
        mask = roipoly(mammo, seg.breast_border(:,1), seg.breast_border(:,2));

        save([mask_dir m_list(ii).name(1:end-4) '_mask.mat'], 'mask');
    end
end
%%
mkdir C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals_half\
m_list = dir('C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\mat\*.mat');
for ii = 1:length(m_list);
    mammo = u_load(['C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\mat\' m_list(ii).name]);
    mammo = imresize(mammo, 0.5, 'bilinear');
    save(['C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals_half\' m_list(ii).name], 'mammo');
    clear mammo;
end
%%
i_list = dir('C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals_half\*.mat');
m_list = dir('C:\isbe\asymmetry_project\data\masks\2004_screening\abnormals\*.mat');
args.num_samples = 1e3;
args.image_dir = 'C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals_half\';
args.mask_dir = 'C:\isbe\asymmetry_project\data\masks\2004_screening\abnormals\';
args.win_size = 3;
args.num_levels = 3;
args.feature_type = 'all';
args.do_max = 0;
args.image_list = i_list(1:10);
args.mask_list = m_list(1:10);
args.save_path = [];
args.plot = 0;
[training_data] = sample_mammo_training_data(args);