%Baby, make me a spiral
bg = u_load('C:\isbe\asymmetry_project\data\synthetic_backgrounds\smooth512\train\bg00001.mat');
label = zeros(512);
label_orientation = zeros(512);
for ii = 0:12:168
    [bar, labeli] = create_sin_bar(4, 4, ii, 512, 512, 0, 256, 256);
    bg = bg + bar;
    label = label | labeli;
    label_orientation(labeli) = ii;
end

figure; imagesc(bg); axis image;
%%
save('C:\isbe\asymmetry_project\data\synthetic_lines\spiral512\spiral001.mat', 'bg');
%%
classify_image_set('191905', 'synthetic_lines\spiral512'); 
%%
classify_image_set('191934', 'synthetic_lines\spiral512', 'forest_dir', 'line_orientation_rfs');
%%
line_prob = u_load('C:\isbe\asymmetry_project\data\synthetic_lines\spiral512\results\191905\spiral001_class.mat');
line_ori = u_load('C:\isbe\asymmetry_project\data\synthetic_lines\spiral512\results\191934\spiral001_class.mat');

for g_width = [8 16 32 64 128]
    [angle_bands dist_sum] = radial_line_projection(...
        line_prob, line_ori, [36 1],...
        fspecial('gaussian', [1 5*g_width], g_width));
    
    figure; imagesc(dist_sum); axis image; colormap(jet(256));
end
%%
g_width = 128;
[angle_bands] = radial_line_projection(...
    line_prob, line_ori, [36 12],...
    fspecial('gaussian', [1 5*g_width], g_width));

%%
figure; imagesc(sum(angle_bands, 3)); axis image; colormap(jet(256));

for ii = 1:12
    figure; imagesc(angle_bands(:,:,ii)); axis image; colormap(jet(256));
end
%%
g_width = 128;
[angle_hist] = radial_line_histogram(...
    line_prob, line_ori, [36 12],...
    fspecial('gaussian', [1 5*g_width], g_width));
%%
for ii = 1:12
    figure; imagesc(angle_hist(:,:,ii)); axis image; colormap(jet(256));
end
%%
mammo = u_load('C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\o04_024RML.mat');
line_map = u_load('C:\isbe\asymmetry_project\data\line_maps\2004_screening\abnormals\024RML_data.mat');
ori_map = u_load('C:\isbe\asymmetry_project\data\orientation_maps\2004_screening\abnormals\024RML_data.mat');

mammo = mammo(1650:1650+511, 1400:1400+511);
line_map = line_map(1650:1650+511, 1400:1400+511);
ori_map = ori_map(1650:1650+511, 1400:1400+511);
pack;

figure; imagesc(mammo); axis image;
figure; imagesc(line_map); axis image;
figure; imagesc(ori_map); axis image; colormap(jet(256));

g_width = 128;
[angle_bands] = radial_line_projection(...
        line_map, ori_map, [36 12],...
        fspecial('gaussian', [1 5*g_width], g_width));
%%
[angle_hist] = radial_line_histogram(...
    line_map, ori_map, [36 12],...
    fspecial('gaussian', [1 5*g_width], g_width));
%%
for ii = 1:12
    figure; imagesc(angle_bands(:,:,ii)); axis image; colormap(jet(256));
    title(['Angle band, sub-band' num2str(ii)]);
%     figure; imagesc(angle_hist(:,:,ii)); axis image; colormap(jet(256));
%     title(['Angle hist, sub-band' num2str(ii)]);
end
figure; imagesc(sum(angle_bands, 3)); axis image; colormap(jet(256));
%%
for ii = 1:12
    figure; imagesc(imfilter(angle_bands(:,:,ii), fspecial('disk', 8))); axis image; colormap(jet(256));
    title(['Angle band, sub-band' num2str(ii)]);
%     figure; imagesc(angle_hist(:,:,ii)); axis image; colormap(jet(256));
%     title(['Angle hist, sub-band' num2str(ii)]);
end
figure; imagesc(imfilter(sum(angle_bands, 3), fspecial('disk', 8))); axis image; colormap(jet(256));