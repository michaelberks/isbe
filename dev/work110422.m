for ii = 5:10:180
    test_bar = create_sin_bar(4, 1, ii, 256, 256, 0.5, 128, 128);
    [line_map, ori_map, scale_map] = gaussian_clover_line(test_bar, [1 2 4 8]);
    figure; imagesc(abs(line_map)); axis image; colormap(gray(256)); hold on;
    quiver(1:4:256, 1:4:256, cos(ori_map(1:4:256,1:4:256)), -sin(ori_map(1:4:256,1:4:256)));
end
%%
for ii = 1:8
    test_bar = create_gauss_bar(ii/2, 1, 23, 256, 256, 128, 128);
    
    figure; 
    for jj = 1:4
        [line_map, ori_map] = gaussian_clover_line(test_bar, 2^(jj-1));
        subplot(2,2,jj);
        imagesc(abs(line_map)); axis image; colormap(gray(256)); hold on; colorbar;
        quiver(1:4:256, 1:4:256, cos(ori_map(1:4:256,1:4:256)), -sin(ori_map(1:4:256,1:4:256)));      
    end
end
%%
for ii = 2:16
    test_bar = create_gauss_bar(ii/2, 1, 23, 256, 256, 128, 128);
    [line_map, ori_map, scale_map] = gaussian_clover_line(test_bar, [1 2 4 8]);
    figure; 
    subplot(1,2,1);
    imagesc(abs(line_map)); axis image; colormap(gray(256)); hold on; colorbar;
    quiver(1:4:256, 1:4:256, cos(ori_map(1:4:256,1:4:256)), -sin(ori_map(1:4:256,1:4:256)));
    subplot(1,2,2);
    imagesc(abs(scale_map)); axis image; colormap(gray(256)); hold on; colorbar;
    
    display(['Mag = ' num2str(abs(line_map(128,128))) ' Scale = ' num2str(scale_map(128,128))]);
end
%%
mkdir Z:\data\synthetic_lines\real512\results\clover\orientations\
for ii = 1:100
    load(['Z:\data\synthetic_lines\real512\image' zerostr(ii,3) '.mat']);
    [line_map, ori_map] = gaussian_clover_line(test_image, [1 2 4 8]);
    save_uint8(['Z:\data\synthetic_lines\real512\results\clover\orientations\image' zerostr(ii,3) '_class.mat'], ori_map);
end
%%
for ii = 1:20
    ret = u_load(['Z:\data\retinograms\DRIVE\test\images_extended\' zerostr(ii,2) '_test_ext.mat']);
    ret = rgb2gray(ret);
    [r c] = size(ret);
    [line_map, ori_map, scale_map] = gaussian_clover_line(ret, [1 2 4 8]);
    ori_complex = abs(line_map).*exp(i*ori_map);
    
    save_uint8(['Z:\data\retinograms\DRIVE\test\orientations\clover\' zerostr(ii,3) '_test_ori.mat'], ori_map);
    
    figure; 
    subplot(1,2,1);
    image(complex2rgb(ori_complex.^2)); axis image;
    subplot(1,2,2);
    imagesc(log(scale_map)); axis image; colormap([1 0 0; 0 1 0; 0 0 1; 0 0 0]);
end
%%
test_dir = 'C:\isbe\asymmetry_project\data\synthetic_lines\real512\';
prob_dir = 'C:\isbe\asymmetry_project\data\synthetic_lines\real512\results\';

wave = 4;
level = 4;

param_dir = ['monogenic_0p65_' zerostr(wave,2) '_' zerostr(level,2) '_ori'];
mkdir([prob_dir param_dir]);

for ii = 1:100
    display(['Processing image ' num2str(ii)]);
    load([test_dir '\image' zerostr(ii,3) '.mat']);

    %Levels 1, min wave length 1, Onf 0.65
    [d d orientation_image] = monogenic_multiscale(test_image, level, wave, 2, 4.65);
    save([prob_dir param_dir '\image' zerostr(ii,3) '_ori.mat'], 'orientation_image');

end           

%%
mkdir Z:\data\retinograms\DRIVE\test\orientations\mono\
for ii = 1:20
    ret = u_load(['Z:\data\retinograms\DRIVE\test\images_extended\' zerostr(ii,2) '_test_ext.mat']);
    ret = rgb2gray(ret);
    [r c] = size(ret);
    [line_map, d, ori_map] = monogenic_multiscale(ret, level, wave, 2, 4.65);
    ori_complex = abs(line_map).*exp(i*ori_map);
    
    save_uint8(['Z:\data\retinograms\DRIVE\test\orientations\mono\' zerostr(ii,3) '_test_ori.mat'], ori_map);
    
    figure; 
    image(complex2rgb(ori_complex.^2)); axis image;

end
%%
mkdir Z:\asym\data\retinograms\DRIVE\test\images_extended\results\g1d\
for ii = 1:20
    ret = u_load(['Z:\asym\data\retinograms\DRIVE\test\images_extended\' zerostr(ii,2) '_test_ext.mat']);
    ret = rgb2gray(ret);
    [line_map, ori_map] = gaussian_1st_derivative_gradient2(ret, [1 2 4 8]);
    ori_complex = abs(line_map).*exp(i*ori_map);
    
    save_uint8(['Z:\asym\data\retinograms\DRIVE\test\images_extended\results\g1d\' zerostr(ii,3) '_test_ori.mat'], ori_map);
    
    figure; 
    image(complex2rgb(ori_complex.^2)); axis image;

end

 
