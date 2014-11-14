for g_width = [8 16 32 64]

    for f_width = [1 4 8 16]

        map_l = u_load(['C:\isbe\dev\ad\mammograms\024LML_dist_map_' zerostr(g_width,3) '_' zerostr(f_width,3) '.mat']);
        map_r = u_load(['C:\isbe\dev\ad\mammograms\024RML_dist_map_' zerostr(g_width,3) '_' zerostr(f_width,3) '.mat']);
        clims = [min([min(map_l(:)) min(map_r(:))]) max([max(map_l(:)) max(map_r(:))])];
        write_im_from_colormap(map_l, ['C:\isbe\dev\ad\mammograms\024LML_dist_map_' zerostr(g_width,3) '_' zerostr(f_width,3) '.jpg'], jet(256), clims);
        write_im_from_colormap(map_r, ['C:\isbe\dev\ad\mammograms\024RML_dist_map_' zerostr(g_width,3) '_' zerostr(f_width,3) '.jpg'], jet(256), clims);
    end
end
%%
i1 = double(imresize(imread('C:\isbe\mammograms\new_CAD\bMP_2004\o04_024LML.bmp'), [1024 NaN], 'bilinear'));
i2 = double(imresize(imread('C:\isbe\mammograms\new_CAD\bMP_2004\o04_024RML.bmp'), [1024 NaN], 'bilinear'));
write_im_from_colormap(i1, 'M:\asymmetry_project\results\contralateral_data\rfs\figures\024RML_cancer.jpg', jet(256));
write_im_from_colormap(i2, 'M:\asymmetry_project\results\contralateral_data\rfs\figures\024LML_cancer.jpg', jet(256));
%%
im1_votes1 = random_forest.image1_votes1 ./ random_forest.image1_total_votes;
im1_votes1(~random_forest.image1_total_votes) = 0;


write_im_from_colormap(im1_votes1, 'M:\asymmetry_project\results\contralateral_data\rfs\figures\024LML_cancer_prob_map.jpg', jet(256));

f1 = figure; hist(im1_votes1(random_forest.image1_total_votes>0), 100);
saveas(f1, 'M:\asymmetry_project\results\contralateral_data\rfs\figures\024LML_cancer_prob_hist.jpg');
%%
del_list = dir('\\isbe-san1\mberks\dev\ad\directed_sums\*008.mat');
for ii = 1:length(del_list)
    delete(['\\isbe-san1\mberks\dev\ad\directed_sums\' del_list(ii).name]);
end
%%
con_list = dir('M:\asymmetry_project\data\mammograms\2004_screening\contralateral_abnormal_roi\*.mat');
for ii = 1:length(con_list)
    a_name = con_list(ii).name(1:6);
    n_name = a_name;
    if strcmpi(a_name(4), 'L')
        n_name(4) = 'R';
    else
        n_name(4) = 'L';
    end
    a_name = dir(['C:\isbe\mammograms\new_CAD\BMP_2004_half\*' a_name '*.mat']);
    n_name = dir(['C:\isbe\mammograms\new_CAD\BMP_2004_half\*' n_name '*.mat']);
    
    a_name_c = ['C:\isbe\mammograms\new_CAD\BMP_2004_half\' a_name(1).name];
    n_name_c = ['C:\isbe\mammograms\new_CAD\BMP_2004_half\' n_name(1).name];
    
    a_name_m = ['M:\asymmetry_project\data\mammograms\2004_screening\abnormals\mat\' a_name(1).name];
    n_name_m = ['M:\asymmetry_project\data\mammograms\2004_screening\abnormals\mat\' n_name(1).name];
    
    copyfile(a_name_c, a_name_m);
    copyfile(n_name_c, n_name_m);
end
%%
con_list = dir('M:\asymmetry_project\data\mammograms\2004_screening\abnormals\mat\*.mat');
for ii = 1:length(con_list)
    m_name = con_list(ii).name(5:10);
    d_names = dir(['\\isbe-san1\mberks\dev\ad\directed_sums\*' m_name '*.mat']);
    for jj = 1:length(d_names)
        old_name = ['\\isbe-san1\mberks\dev\ad\directed_sums\' d_names(jj).name];
        new_name = ['M:\asymmetry_project\data\radial_maps\2004_screening\abnormals\' d_names(jj).name(1:7) 'rad_map' d_names(jj).name(15:18) '.mat'];
        
        if ~exist(new_name, 'file')
            copyfile(old_name, new_name);
        end
    end
end
%%
norm_list = dir('C:\isbe\mammograms\new_CAD\BMP_2004_Normals\*.bmp');
mkdir C:\isbe\mammograms\new_CAD\BMP_2004_Normals_half;
mkdir M:\asymmetry_project\data\mammograms\2004_screening\normals\mat\;
for ii = 1:length(norm_list)
    mammo = imread(['C:\isbe\mammograms\new_CAD\BMP_2004_Normals\' norm_list(ii).name]);
    mammo_small = imresize(mammo, 0.5, 'bilinear'); clear mammo;
    
    save(['C:\isbe\mammograms\new_CAD\BMP_2004_Normals_half\' norm_list(ii).name(1:end-3) 'mat'], 'mammo_small');
    save(['M:\asymmetry_project\data\mammograms\2004_screening\normals\mat\' norm_list(ii).name(1:end-3) 'mat'], 'mammo_small');
    clear mammo_small;
end
%%
con_list = dir('M:\asymmetry_project\data\mammograms\2004_screening\abnormals\mat\*.mat');
missing = [];
jj = 1;
for ii = 1:length(con_list)
    m_name = con_list(ii).name(5:10);
    d_names = dir(['M:\asymmetry_project\data\radial_maps\2004_screening\abnormals\*' m_name '*.mat']);
    if length(d_names) ~= 4
        display([m_name ' has ' num2str(length(d_names)) ' maps']);
        missing{jj,1} = m_name;
        jj = jj+1;
    end
end
%%
for ii = 1:size(missing,1);
    try
        load(['\\Isbe-san1\mberks\dev\ad\' missing{ii} '_data.mat'], 'line_ori', 'line_prob');
    catch
        display(['Couldn''t find \\Isbe-san1\mberks\dev\ad\' missing{ii} '_data.mat']);
        continue;
    end
    line_ori = imresize(line_ori, 0.5, 'bilinear');
    line_prob = imresize(line_prob, 0.5, 'bilinear');
    
    for g_width = [8 16 32 64]
        fname = ['M:\asymmetry_project\data\radial_maps\2004_screening\abnormals\' missing{ii} '_rad_map_' zerostr(g_width,3) '_001.mat'];
        [angle_bands dist_sum] = radial_line_projection(line_prob, line_ori, [36 1], fspecial('gaussian', [1 5*g_width], g_width)); %#ok
        save(fname, 'dist_sum');
    end
end
%%
con_list = dir('M:\asymmetry_project\data\mammograms\2004_screening\abnormals\mat\*.mat');
% mkdir C:\isbe\asymmetry_project\data\line_detection_maps\2004_screening\abnormals\
% mkdir C:\isbe\asymmetry_project\data\orientation_maps\2004_screening\abnormals\

for ii = 1:length(con_list)
    m_name = con_list(ii).name(5:10);
    l_name = ['C:\isbe\asymmetry_project\data\line_detection_maps\2004_screening\abnormals\' m_name '_data.mat'];
    o_name = ['C:\isbe\asymmetry_project\data\orientation_maps\2004_screening\abnormals\' m_name '_data.mat'];
    
    if ~exist(l_name, 'file')
        try
            load(['C:\isbe\asymmetry_project\misc\' m_name '_data.mat'], 'line_ori', 'line_prob');
        catch
            display(['Couldn''t find C:\isbe\asymmetry_project\misc\' m_name '_data.mat']);
            continue;
        end
        save(l_name, 'line_prob');
        save(o_name, 'line_ori');
        delete(['C:\isbe\asymmetry_project\misc\' m_name '_data.mat']);
        clear line_*
    end    
end
%%