%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Script describing how various mammographic datasets have been created for
% the asymmetry project
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%%
%--------------------------------------------------------------------------
% 1: Temporal dataset
%--------------------------------------------------------------------------
abnormal_list = dir('I:\mammograms\abnormal cases\*.tif');

for ii = 0 + (1:20)
    mam = imread(['I:\mammograms\abnormal cases\' abnormal_list(ii).name]);
    figure; imagesc(mam); axis image; colormap(gray(256));
    title(abnormal_list(ii).name);
    clear mam;
end
%%
abnormal_list = dir('I:\mammograms\abnormal cases\*.tif');
mkdir M:\asymmetry_project\data\mammograms\marco\abnormals
%%
for ii = 1:length(abnormal_list);
    mam = imread(['I:\mammograms\abnormal cases\' abnormal_list(ii).name]);
    mam = uint8(mam / 256);
    mam_small = imresize(mam, 0.5, 'bilinear');
    clear mam;
    
    imwrite(mam_small, ['M:\asymmetry_project\data\mammograms\marco\abnormals\' abnormal_list(ii).name(1:end-3) 'jpg']);
    clear mam_small;
end
%%
normal_list = dir('I:\mammograms\normal cases\*.tif');
mkdir M:\asymmetry_project\data\mammograms\marco\normals
%%
for ii = 1:length(normal_list);
    mam = imread(['I:\mammograms\normal cases\' normal_list(ii).name]);
    mam_small = imresize(mam, 0.5, 'bilinear');
    clear mam;
    
    imwrite(mam_small, ['M:\asymmetry_project\data\mammograms\marco\normals\' normal_list(ii).name(1:end-3) 'jpg'], 'bitdepth', 12);
    clear mam_small;
end
%%
%--------------------------------------------------------------------------
% 2: MB's screening cancers dataset
%--------------------------------------------------------------------------

% First check the list of masses and discard cases with a mass in both
% breasts
mass_list = dir('C:\isbe\dev\annotations\*.mat');
pairs = false(length(mass_list),1);
for ii = 1:length(mass_list)
    
    mass_name = mass_list(ii).name;
    opp_name = mass_name;
    if ~isempty(strfind(mass_name, 'R'))
        opp_name(9) = 'L';
    else
        opp_name(9) = 'R';
    end
    
    for jj = 1:length(mass_list)
        if strcmp(opp_name, mass_list(jj).name)
            display([mass_name ' and ' opp_name ' are a pair']);
            pairs(ii) = 1;
            break;
        end
    end
end
mass_list(pairs) = [];
%%
failed_list = [];
% Now for each mass, extract
for ii = 1:length(mass_list)
    
    try
        mass = u_load(['C:\isbe\dev\annotations\' mass_list(ii).name]);
        mass_name = mass.name;
        opp_name = mass_name;
        if strcmpi(mass_name(end-6), 'R')
            opp_name(end-6) = 'L';
        else
            opp_name(end-6) = 'R';
        end

        mammo1 = imread(mass_name);
        mammo2 = imread(opp_name);

        seg1 = u_load(['C:\isbe\dev\segmentation\breast_borders\' mass_name(37:46) '_segmentation.mat']);
        seg2 = u_load(['C:\isbe\dev\segmentation\breast_borders\' opp_name(37:46) '_segmentation.mat']);

        [seg1.breast_border seg1.breast_air] = segment_breast_resize(size(mammo1), seg1);
        [seg2.breast_border seg2.breast_air] = segment_breast_resize(size(mammo2), seg2);

        source_breast = seg1.breast_border(seg1.breast_air,:);
        target_breast = seg2.breast_border(seg2.breast_air,:);

        source_pos = [mass.C2 mass.R1];
        target_pos = round(select_corresponding_position(source_breast, target_breast, source_pos, 1));

        target_pos(2,1) = target_pos(1,1) + mass.C2 - mass.C1;
        target_pos(2,2) = target_pos(1,2) + mass.R2 - mass.R1;

        target_pos(1,1) = max(target_pos(1,1), 1);
        target_pos(1,2) = max(target_pos(1,2), 1);
        target_pos(2,1) = min(target_pos(2,1), size(mammo2,2));
        target_pos(2,2) = min(target_pos(2,2), size(mammo2,1));

    %     figure; 
    %     subplot(1,2,1); imagesc(mammo1); axis image; colormap(gray(256)); hold on;
    %     plot(source_breast(:,1), source_breast(:,2));
    %     plot([mass.C1 mass.C2 mass.C2 mass.C1 mass.C1],...
    %          [mass.R1 mass.R1 mass.R2 mass.R2 mass.R1], 'r');
    %     
    %     subplot(1,2,2); imagesc(mammo2); axis image; colormap(gray(256)); hold on;
    %     plot(target_breast(:,1), target_breast(:,2));
    %     plot([target_pos(1,1) target_pos(2,1) target_pos(2,1) target_pos(1,1) target_pos(1,1)],...
    %          [target_pos(1,2) target_pos(1,2) target_pos(2,2) target_pos(2,2) target_pos(1,2)], 'g');
    %     
        contralateral_pair.abnormal_roi = mammo1(mass.R1:mass.R2, mass.C1:mass.C2);
        contralateral_pair.normal_roi = mammo2(target_pos(1,2):target_pos(2,2), target_pos(1,1):target_pos(2,1));
        clear mammo* seg*

        contralateral_pair.right = strcmpi(mass_name(end-6), 'R');
        if contralateral_pair.right
            contralateral_pair.abnormal_roi = fliplr(contralateral_pair.abnormal_roi);
        else
            contralateral_pair.normal_roi = fliplr(contralateral_pair.normal_roi);
        end

        contralateral_pair.abnormal_pos = [mass.C1 mass.R1; mass.C2 mass.R2];
        contralateral_pair.normal_pos = target_pos;
        contralateral_pair.abnormal_name = mass_name;
        contralateral_pair.normal_name = opp_name;

        save(['M:\asymmetry_project\data\contralateral\' mass_list(ii).name(6:end-4) '.mat'], 'contralateral_pair');
        clear contralateral_pair
    catch err
        %add failed mammogram name to failures list
        failed_list(end+1,1).name = mass_list(ii).name; %#ok;
        display(['Skipping ', mass_list(ii).name '. ' err.message]);
    end
%     figure; 
%     subplot(1,2,1); imagesc(abnormal_roi); axis image; colormap(gray(256)); hold on;
%     subplot(1,2,2); imagesc(normal_roi); axis image; colormap(gray(256)); hold on;
    
    
end
%%
% Now do the same again but for the normals

%First segment the images
[failed_list] = segment_breast_batch(...
    'C:\isbe\mammograms\new_CAD\BMP_2004_Normals\',...
    'C:\isbe\mammograms\new_CAD\BMP_2004_Normals\segmentations\',...
    'bmp'); %#ok redo any failures

%%
con_list = dir('M:\asymmetry_project\data\mammograms\2004_screening\contralateral_normal_19rois\*.mat');
for ii = 1:length(con_list)
    load(['M:\asymmetry_project\data\mammograms\2004_screening\contralateral_normal_19rois\' con_list(ii).name]);

% %     contralateral_pair.abnormal_pos = round(contralateral_pair.abnormal_pos);
% %     
% %     if all(diff(contralateral_pair.abnormal_pos)+1 == size(contralateral_pair.abnormal_roi))
% %         contralateral_pair.abnormal_pos = fliplr(contralateral_pair.abnormal_pos);
% %         display(['Flipped ', num2str(ii)]);
% %     elseif all(diff(fliplr(contralateral_pair.abnormal_pos))+1 == size(contralateral_pair.abnormal_roi))
% %         display(['Ok: ', num2str(ii)]);
% %     else
% %         display(['???? ', num2str(ii)]);
% %     end
% %     save(['M:\asymmetry_project\data\mammograms\2004_screening\contralateral_normal_19rois\' con_list(ii).name], 'contralateral_pair');    
% % end
% % %%

    l_name = con_list(ii).name(1:6);
    r_name = l_name;
    r_name(4) = 'R';
    mammo_l_name = ['M:\asymmetry_project\data\mammograms\2004_screening\normals\bmp\' l_name '.bmp'];
    mammo_r_name = ['M:\asymmetry_project\data\mammograms\2004_screening\normals\bmp\' r_name '.bmp'];
    
    mammo_l = imread(mammo_l_name);
    mammo_r = imread(mammo_r_name);

    seg_l = u_load(['M:\asymmetry_project\data\segmentations\2004_screening\normals\' l_name '_segmentation.mat']);
    seg_r = u_load(['M:\asymmetry_project\data\segmentations\2004_screening\normals\' r_name '_segmentation.mat']);
    
    [seg_l.breast_border seg_l.breast_air] = segment_breast_resize(size(mammo_l), seg_l);
    [seg_r.breast_border seg_r.breast_air] = segment_breast_resize(size(mammo_r), seg_r);
    
    source_breast = seg_l.breast_border(seg_l.breast_air,:);
    target_breast = seg_r.breast_border(seg_r.breast_air,:);

    source_pos = [contralateral_pair.abnormal_pos(2,1) contralateral_pair.abnormal_pos(1,2)];
    target_pos = round(select_corresponding_position(source_breast, target_breast, source_pos, 1));

    target_pos(2,1) = target_pos(1,1) + size(contralateral_pair.abnormal_roi,2) - 1;
    target_pos(2,2) = target_pos(1,2) + size(contralateral_pair.abnormal_roi,1) - 1;
    
    try
        contralateral_pair.normal_pos = target_pos;
        contralateral_pair.normal_roi = mammo_r(target_pos(1,2):target_pos(2,2), target_pos(1,1):target_pos(2,1));
        contralateral_pair.normal_roi = fliplr(contralateral_pair.normal_roi);
        contralateral_pair.abnormal_name = mammo_l_name;
        contralateral_pair.normal_name = mammo_r_name;

        save(['M:\asymmetry_project\data\mammograms\2004_screening\contralateral_normal_19rois\' con_list(ii).name], 'contralateral_pair');
    catch
        display(num2str(ii));
        display(num2str(target_pos));
    end
    
%     figure; 
%     subplot(1,2,1); imagesc(mammo_l); axis image; colormap(gray(256)); hold on;
%     plot(seg_l.breast_border(:,1), seg_l.breast_border(:,2), 'b');
%     plot([contralateral_pair.abnormal_pos(1,1) contralateral_pair.abnormal_pos(2,1) contralateral_pair.abnormal_pos(2,1) contralateral_pair.abnormal_pos(1,1) contralateral_pair.abnormal_pos(1,1)],...
%          [contralateral_pair.abnormal_pos(1,2) contralateral_pair.abnormal_pos(1,2) contralateral_pair.abnormal_pos(2,2) contralateral_pair.abnormal_pos(2,2) contralateral_pair.abnormal_pos(1,2)], 'g');
%     subplot(1,2,2); imagesc(mammo_r); axis image; colormap(gray(256)); hold on;
%     plot(seg_r.breast_border(:,1), seg_r.breast_border(:,2), 'b');
%     plot([contralateral_pair.normal_pos(1,1) contralateral_pair.normal_pos(2,1) contralateral_pair.normal_pos(2,1) contralateral_pair.normal_pos(1,1) contralateral_pair.normal_pos(1,1)],...
%          [contralateral_pair.normal_pos(1,2) contralateral_pair.normal_pos(1,2) contralateral_pair.normal_pos(2,2) contralateral_pair.normal_pos(2,2) contralateral_pair.normal_pos(1,2)], 'g');
%     figure; 
%     subplot(1,2,1); imagesc(contralateral_pair.abnormal_roi); axis image; colormap(gray(256)); hold on;
%     subplot(1,2,2); imagesc(contralateral_pair.normal_roi); axis image; colormap(gray(256)); hold on;

    
    


    clear contralateral_pair
    
end
%%
con_list = dir('M:\asymmetry_project\data\mammograms\2004_screening\contralateral_normal_rois\*.mat');
for ii = 1:10
    load(['M:\asymmetry_project\data\mammograms\2004_screening\contralateral_normal_rois\' con_list(ii).name]);

    l_name = con_list(ii).name(1:6);
    r_name = l_name;
    r_name(4) = 'R';
    mammo_l_name = ['M:\asymmetry_project\data\mammograms\2004_screening\normals\bmp\' l_name '.bmp'];
    mammo_r_name = ['M:\asymmetry_project\data\mammograms\2004_screening\normals\bmp\' r_name '.bmp'];
   
    contralateral_pair.abnormal_name = mammo_l_name;
    contralateral_pair.normal_name = mammo_r_name;
    
    save(['M:\asymmetry_project\data\mammograms\2004_screening\contralateral_normal_rois\' con_list(ii).name], 'contralateral_pair');

    clear contralateral_pair
    
end
%%
con_list = dir('C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\*.mat');
for ii = 1:length(con_list)
    m_name = con_list(ii).name(1:6);
    a_name = dir(['D:\isbe\dev\annotations\*' m_name '*.mat']);
    
    if ~isempty(a_name)
        mam = u_load(['C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\' con_list(ii).name]);
        [r c] = size(mam);
        meta_xy = [];
        for jj = 1:length(a_name)
            anno = u_load(['D:\isbe\dev\annotations\' a_name(jj).name]);
            x = (anno.mass_outline(:,1) + anno.C1) / (2*c);
            y = (anno.mass_outline(:,2) + anno.R1) / (2*r);
            meta_xy = [meta_xy; x y]; %#ok
        end
        save(['C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\meta\' m_name '_meta.mat'], 'meta_xy');
    end
end
%%
con_list = dir('C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\*.mat');
mkdir C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\meta2;

for ii = 1:length(con_list)
    m_name = con_list(ii).name(1:6);
    a_name = dir(['D:\isbe\dev\annotations\*' m_name '*.mat']);
    if ~isempty(a_name)
        mam = u_load(['C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\' con_list(ii).name]);
        [r c] = size(mam);
        meta_xy = cell(length(a_name),1);
        
        %figure; hold on;
        for jj = 1:length(a_name)
            anno = u_load(['D:\isbe\dev\annotations\' a_name(jj).name]);
            
            x = (anno.mass_outline(:,1) + anno.C1);
            y = (anno.mass_outline(:,2) + anno.R1);
            
            %plot(x, y, 'b');
            
            x = (x - mean(x))*1.1 + mean(x);
            y = (y - mean(y))*1.1 + mean(y);
            
            
            K = convhull(x,y);
            meta_xy{jj} = [x(K)/(2*c) y(K)/(2*r)];
            
            %plot(meta_xy{jj}(:,1), meta_xy{jj}(:,2), 'r');
            %plot(x, y, 'kx');
        end
        save(['C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\mat\meta2\' m_name '_meta.mat'], 'meta_xy');
    end
end
%%
%--------------------------------------------------------------------------
seg_list = dir('C:\isbe\asymmetry_project\data\segmentations\2004_screening\abnormals\*mat');
mam_list = dir('C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\mat\*mat');

mam_sizes = zeros(292,1);
for ii = 1:292
    mam = u_load(['C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\mat\' mam_list(ii).name]);
    seg = u_load(['C:\isbe\asymmetry_project\data\segmentations\2004_screening\abnormals\' seg_list(ii).name]);
    
    seg.breast_border = segment_breast_resize(size(mam), seg); clear mam;
    
    mam_sizes(ii) = polyarea(seg.breast_border(:,1), seg.breast_border(:,2));
end
%%
con_list = dir('C:\isbe\asymmetry_project\data\mammograms\2004_screening\contralateral_abnormal_roi\*.mat');
roi_sizes = zeros(length(con_list), 1); 
for ii = 1:length(con_list)
    con = u_load(['C:\isbe\asymmetry_project\data\mammograms\2004_screening\contralateral_abnormal_roi\' con_list(ii).name]);
    roi_sizes(ii) = numel(con.abnormal_roi);
    clear con;
end
%% Convert segmentations into binary masks
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
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%% Extract a set of real 512x512 regions from the normal mammograms
mb_sample_valid_region(...
    'WindowSize', 512,...
    'MammoDir', 'C:\isbe\asymmetry_project\data\mammograms\2004_screening\normals',...
    'SaveDir', 'C:\isbe\asymmetry_project\data\synthetic_backgrounds\real512',...
    'NumRegions', 5,...
    'NameStem', 'bg',...
    'MammoFiles', [], ...
    'MaskDir', 'C:\isbe\asymmetry_project\data\masks\2004_screening\normals',...
    'MaskFiles',[],...
    'ResizeFactor', 4, ...
    'ImageFormat', '.mat');
%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%% Create ROI of the masses and their associated maps
mass_list = dir('C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\meta\*.mat');
mass_names = get_mammo_info(mass_list);
roi_dirs = {...
    'C:\isbe\asymmetry_project\data\mammograms\',...
    'C:\isbe\asymmetry_project\data\line_maps\rf\',...
    'C:\isbe\asymmetry_project\data\orientation_maps\rf\',...
    'C:\isbe\asymmetry_project\data\orientation_maps\g2d\',...
    'C:\isbe\asymmetry_project\data\line_maps\g2d\',...
    'Z:\data\orientation_maps\rf\',...
    'Z:\data\orientation_maps\rf_thin\',...
    'Z:\data\orientation_maps\g2d\',...
    'Z:\data\line_maps\g2d\',...
    'Z:\data\line_maps\rf_prob\',...
    'Z:\data\orientation_maps\rf_prob\',...
    'C:\isbe\asymmetry_project\data\masks\'};
mam_types = {...
    '2004_screening_processed\',...
    '2004_screening\',...
    '2004_screening_processed\',...
    '2004_screening_processed\',...
    '2004_screening_processed\',...
    '2004_screening_processed\',...
    '2004_screening_processed\',...
    '2004_screening_processed\',...
    '2004_screening_processed\',...
    '2004_screening_processed\',...
    '2004_screening_processed\',...
    '2004_screening_processed\'};
win_sizes = 800*ones(length(roi_dirs),1);
for ii = 12:12%2:length(roi_dirs);
    roi_names = match_mammo_names([roi_dirs{ii} mam_types{ii} 'abnormals\'], mass_names);
    mkdir([roi_dirs{ii} mam_types{ii} 'mass_roi\']);
    if ii == 1
        mkdir([roi_dirs{ii} mam_types{ii} 'mass_roi\meta\']);
    end
    for jj = 1:length(mass_names);
        mass_xy = u_load(['C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\meta\' mass_list(jj).name]);
        full_image = load_uint8([roi_dirs{ii} mam_types{ii} 'abnormals\' roi_names{jj}]);
        [r c] = size(full_image);
        x = mass_xy(:,1) * c;
        y = mass_xy(:,2) * r;
        xc = round(mean(x));
        yc = round(mean(y));
        roi = sample_window(full_image, win_sizes(ii), yc, xc,0);
        
        save_uint8([roi_dirs{ii} mam_types{ii} 'mass_roi\' mass_names{jj} '_roi.mat'], roi);
        if ii == 1
            mass_xy(:,1) = (x - xc + win_sizes(ii)/2) / win_sizes(ii);
            mass_xy(:,2) = (y - yc + win_sizes(ii)/2) / win_sizes(ii);
            save([roi_dirs{ii} mam_types{ii} 'mass_roi\meta\' mass_names{jj} '_meta.mat'], 'mass_xy');
        end
        if jj < 4
            figure; imagesc(complex2rgb(roi)); axis image;
        end
    end
end
%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%Create mass BGs

mkdir C:\isbe\asymmetry_project\data\synthetic_backgrounds\mass_roi\train\
mkdir C:\isbe\asymmetry_project\data\synthetic_backgrounds\mass_roi\masks\
for ii = 1:146
    
    %Copy mass roi into mass dir
    copyfile(...
        ['C:\isbe\asymmetry_project\data\mammograms\2004_screening_processed\mass_roi\' mass_names{ii} '_roi.mat'],...
        ['C:\isbe\asymmetry_project\data\synthetic_backgrounds\mass_roi\train\bg' zerostr(ii,5) '.mat']);
    %Copy mass mask into mask dir
    copyfile(...
        ['C:\isbe\asymmetry_project\data\masks\2004_screening_processed\mass_roi\' mass_names{ii} '_roi.mat'],...
        ['C:\isbe\asymmetry_project\data\synthetic_backgrounds\mass_roi\masks\bg' zerostr(ii,5) '_mask.mat']);
    %Copy existing bg into mass dir
    copyfile(...
        ['C:\isbe\asymmetry_project\data\synthetic_backgrounds\real512\train\bg' zerostr(ii,5) '.mat'],...
        ['C:\isbe\asymmetry_project\data\synthetic_backgrounds\mass_roi\train\bg' zerostr(ii+146,5) '.mat']);
end
        
