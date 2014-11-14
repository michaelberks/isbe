%mkdir C:\isbe\asymmetry_project\data\mammograms\2004_screening\contralateral_abnormal_roi2\;

mkdir C:\isbe\asymmetry_project\data\mammograms\2004_screening\contralateral_roi\abnormals\
mkdir C:\isbe\asymmetry_project\data\mammograms\2004_screening\contralateral_roi\normals\
mkdir C:\isbe\asymmetry_project\data\masks\2004_screening\contralateral_roi\abnormals\
mkdir C:\isbe\asymmetry_project\data\masks\2004_screening\contralateral_roi\normals\

roi_list = dir('C:\isbe\asymmetry_project\data\mammograms\2004_screening\contralateral_abnormal_roi\*.mat');
ab_names = get_mammo_info(roi_list);
for ii = 1:length(roi_list);
    load(['C:\isbe\asymmetry_project\data\mammograms\2004_screening\contralateral_abnormal_roi\' roi_list(ii).name]);
    ab_mammo = u_load(['C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\mat\o04_' ab_names{ii} '.mat']);
    ab_mask = u_load(['C:\isbe\asymmetry_project\data\masks\2004_screening\abnormals\o04_' ab_names{ii} '_mask.mat']);
    
    contralateral_pair.abnormal_pos = round(contralateral_pair.abnormal_pos/2);
    
    ra1 = contralateral_pair.abnormal_pos(1,2);
    ra2 = contralateral_pair.abnormal_pos(2,2);
    ca1 = contralateral_pair.abnormal_pos(1,1);
    ca2 = contralateral_pair.abnormal_pos(2,1);
    abnormal_roi = ab_mammo(ra1:ra2, ca1:ca2);
    abnormal_mask = ab_mask(ra1:ra2, ca1:ca2);
    
    norm_name = ab_names{ii};
    if contralateral_pair.right
        norm_name(4) = 'L';
    else
        norm_name(4) = 'R';
    end
    
    norm_mammo = u_load(['C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\mat\o04_' norm_name '.mat']);
    norm_mask = u_load(['C:\isbe\asymmetry_project\data\masks\2004_screening\abnormals\o04_' norm_name '_mask.mat']);
    
    contralateral_pair.normal_pos = round(contralateral_pair.normal_pos/2);
    
    rn1 = contralateral_pair.normal_pos(1,2);
    rn2 = contralateral_pair.normal_pos(2,2);
    cn1 = contralateral_pair.normal_pos(1,1);
    cn2 = contralateral_pair.normal_pos(2,1);
    normal_roi = norm_mammo(rn1:rn2, cn1:cn2);
    normal_mask = norm_mask(rn1:rn2, cn1:cn2);
    
%     if contralateral_pair.right
%         abnormal_roi = fliplr(abnormal_roi);
%         abnormal_mask = fliplr(abnormal_mask);
%     else
%         normal_roi = fliplr(normal_roi);
%         normal_mask = fliplr(normal_mask);
%     end
%     
%     figure; 
%     subplot(2,2,1); imagesc(abnormal_roi); axis image;
%     subplot(2,2,2); imagesc(normal_roi); axis image;
%     subplot(2,2,3); imagesc(abnormal_mask); axis image;
%     subplot(2,2,4); imagesc(normal_mask); axis image;
    save(['C:\isbe\asymmetry_project\data\mammograms\2004_screening\contralateral_roi\abnormals\' ab_names{ii} '.mat'], 'abnormal_roi');
    save(['C:\isbe\asymmetry_project\data\mammograms\2004_screening\contralateral_roi\normals\' norm_name '.mat'], 'normal_roi');
    save(['C:\isbe\asymmetry_project\data\masks\2004_screening\contralateral_roi\abnormals\' ab_names{ii} '_mask.mat'], 'abnormal_mask');
    save(['C:\isbe\asymmetry_project\data\masks\2004_screening\contralateral_roi\normals\' norm_name '_mask.mat'], 'normal_mask');
end
%%
%mkdir C:\isbe\asymmetry_project\data\mammograms\2004_screening\contralateral_abnormal_roi2\;

mkdir Z:\asymmetry_project\data\mammograms\2004_screening\contralateral\abnormals\
mkdir Z:\asymmetry_project\data\mammograms\2004_screening\contralateral\normals\
mkdir Z:\asymmetry_project\data\masks\2004_screening\contralateral\abnormals\
mkdir Z:\asymmetry_project\data\masks\2004_screening\contralateral\normals\

ab_list = dir('C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\mat\meta\*.mat');
ab_names = get_mammo_info(ab_list);
for ii = 1:length(ab_names);
    norm_name = ab_names{ii};
    if strcmpi(ab_names{ii}(4), 'R')
        norm_name(4) = 'L';
    else
        norm_name(4) = 'R';
    end
    copyfile(...
        ['Z:\asymmetry_project\data\mammograms\2004_screening\abnormals\mat\o04_' ab_names{ii} '.mat'],...
        ['Z:\asymmetry_project\data\mammograms\2004_screening\contralateral\abnormals\' ab_names{ii} '.mat']);
    copyfile(...
        ['Z:\asymmetry_project\data\mammograms\2004_screening\abnormals\mat\o04_' norm_name '.mat'],...
        ['Z:\asymmetry_project\data\mammograms\2004_screening\contralateral\normals\' norm_name '.mat']);
    
    copyfile(...
        ['Z:\asymmetry_project\data\masks\2004_screening\abnormals\o04_' ab_names{ii} '_mask.mat'],...
        ['Z:\asymmetry_project\data\masks\2004_screening\contralateral\abnormals\' ab_names{ii} '_mask.mat']);
    copyfile(...
        ['Z:\asymmetry_project\data\masks\2004_screening\abnormals\o04_' norm_name '_mask.mat'],...
        ['Z:\asymmetry_project\data\masks\2004_screening\contralateral\normals\' norm_name '_mask.mat']);
end
%%
clear;
mammo = imresize(u_load('Z:\asymmetry_project\data\mammograms\2004_screening\contralateral\abnormals\024RML.mat'), 0.5, 'bilinear');
line_map = imresize(u_load('C:\isbe\asymmetry_project\data\line_maps\2004_screening\abnormals\024RML_data.mat'), 0.5, 'bilinear');
ori_map = imresize(u_load('C:\isbe\asymmetry_project\data\orientation_maps\2004_screening\abnormals\024RML_data.mat'), 0.5, 'bilinear');
abnormal_map = fliplr(imresize(1-u_load('Z:\asymmetry_project\data\mammograms\2004_screening\contralateral\abnormals\results\abnormal_maps\024RML_class.mat'), 0.5, 'bilinear'));
%
figure; imagesc(mammo); axis image;
figure; imagesc(line_map); axis image;
figure; imagesc(line_map .* abnormal_map); axis image;
%
g_width = 64;
[rad_map_inf_1 rad_map128_1] = radial_line_projection(line_map, ori_map, [36 1], fspecial('gaussian', [1 5*g_width], g_width));
[rad_map_inf_2 rad_map128_2] = radial_line_projection(line_map .* abnormal_map, ori_map, [36 1], fspecial('gaussian', [1 5*g_width], g_width));

figure; imagesc(imfilter(rad_map128_1, fspecial('gaussian', 40, 4))); axis image; colormap(jet(256));
figure; imagesc(imfilter(rad_map128_2, fspecial('gaussian', 40, 4))); axis image; colormap(jet(256)); 
%%
clear;
mammo = imresize(u_load('Z:\asymmetry_project\data\mammograms\2004_screening\contralateral\normals\024LML.mat'), 0.5, 'bilinear');
line_map = imresize(u_load('C:\isbe\asymmetry_project\data\line_maps\2004_screening\abnormals\024LML_data.mat'), 0.5, 'bilinear');
ori_map = imresize(u_load('C:\isbe\asymmetry_project\data\orientation_maps\2004_screening\abnormals\024LML_data.mat'), 0.5, 'bilinear');
abnormal_map = fliplr(imresize(1-u_load('Z:\asymmetry_project\data\mammograms\2004_screening\contralateral\normals\results\abnormal_maps\024LML_class.mat'), 0.5, 'bilinear'));
%
figure; imagesc(mammo); axis image;
figure; imagesc(line_map); axis image;
figure; imagesc(line_map .* abnormal_map); axis image;
%
g_width = 64;
[rad_map_inf_1 rad_map128_1] = radial_line_projection(line_map, ori_map, [36 1], fspecial('gaussian', [1 5*g_width], g_width));
[rad_map_inf_2 rad_map128_2] = radial_line_projection(line_map .* abnormal_map, ori_map, [36 1], fspecial('gaussian', [1 5*g_width], g_width));

figure; imagesc(imfilter(rad_map128_1, fspecial('gaussian', 40, 4))); axis image; colormap(jet(256));
figure; imagesc(imfilter(rad_map128_2, fspecial('gaussian', 40, 4))); axis image; colormap(jet(256)); 
%%
mkdir C:\isbe\asymmetry_project\data\synthetic_lines\lines512\labels
mkdir C:\isbe\asymmetry_project\data\synthetic_lines\curves512\labels
for ii = 1:100
    load(['C:\isbe\asymmetry_project\data\synthetic_lines\lines512\image' zerostr(ii,3)]);
    save(['C:\isbe\asymmetry_project\data\synthetic_lines\lines512\labels\label' zerostr(ii,3)], 'label*');
    save(['C:\isbe\asymmetry_project\data\synthetic_lines\lines512\image' zerostr(ii,3)], 'test_image');
    clear label* test_image
    
    load(['C:\isbe\asymmetry_project\data\synthetic_lines\curves512\image' zerostr(ii,3)]);
    save(['C:\isbe\asymmetry_project\data\synthetic_lines\curves512\labels\label' zerostr(ii,3)], 'label*');
    save(['C:\isbe\asymmetry_project\data\synthetic_lines\curves512\image' zerostr(ii,3)], 'test_image');
    clear label* test_image
end
%%
mkdir C:\isbe\asymmetry_project\data\masks\2004_screening\contralateral_roi\abnormal_lines\
mkdir C:\isbe\asymmetry_project\data\masks\2004_screening\contralateral_roi\normal_lines\

roi_list = dir('C:\isbe\asymmetry_project\data\mammograms\2004_screening\contralateral_abnormal_roi\*.mat');
ab_names = get_mammo_info(roi_list);
for ii = 1:length(roi_list);
    load(['C:\isbe\asymmetry_project\data\mammograms\2004_screening\contralateral_abnormal_roi\' roi_list(ii).name]);
    
    ab_line = u_load(['C:\isbe\asymmetry_project\data\line_maps\2004_screening\abnormals\' ab_names{ii} '_data.mat']);
    
    contralateral_pair.abnormal_pos = round(contralateral_pair.abnormal_pos/2);
    
    ra1 = contralateral_pair.abnormal_pos(1,2);
    ra2 = contralateral_pair.abnormal_pos(2,2);
    ca1 = contralateral_pair.abnormal_pos(1,1);
    ca2 = contralateral_pair.abnormal_pos(2,1);

    abnormal_mask = ab_line(ra1:ra2, ca1:ca2) > 0.5;
    
    norm_name = ab_names{ii};
    if contralateral_pair.right
        norm_name(4) = 'L';
    else
        norm_name(4) = 'R';
    end
    
    norm_line = u_load(['C:\isbe\asymmetry_project\data\line_maps\2004_screening\abnormals\' norm_name '_data.mat']);
    
    contralateral_pair.normal_pos = round(contralateral_pair.normal_pos/2);
    
    rn1 = contralateral_pair.normal_pos(1,2);
    rn2 = contralateral_pair.normal_pos(2,2);
    cn1 = contralateral_pair.normal_pos(1,1);
    cn2 = contralateral_pair.normal_pos(2,1);

    normal_mask = norm_line(rn1:rn2, cn1:cn2) > 0.5;
    
%     figure; 
%     subplot(2,2,1); imagesc(contralateral_pair.abnormal_roi); axis image;
%     subplot(2,2,2); imagesc(contralateral_pair.normal_roi); axis image;
%     subplot(2,2,3); imagesc(abnormal_mask); axis image;
%     subplot(2,2,4); imagesc(normal_mask); axis image;
    
    save(['C:\isbe\asymmetry_project\data\masks\2004_screening\contralateral_roi\abnormal_lines\' ab_names{ii} '_mask.mat'], 'abnormal_mask');
    save(['C:\isbe\asymmetry_project\data\masks\2004_screening\contralateral_roi\normal_lines\' norm_name '_mask.mat'], 'normal_mask');
end
