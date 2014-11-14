mammo = imresize(u_load('C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\o04_024RML.mat'), 0.5, 'bilinear');
line_map = imresize(u_load('C:\isbe\asymmetry_project\data\line_maps\2004_screening\abnormals\024RML_data.mat'), 0.5, 'bilinear');
ori_map = imresize(u_load('C:\isbe\asymmetry_project\data\orientation_maps\2004_screening\abnormals\024RML_data.mat'), 0.5, 'bilinear');

label = mb_non_maximal_supp(line_map, ori_map, 1);

w_ori = line_map .* exp(1i * pi * ori_map / 180);
%%
display_orientation(label, w_ori, 1, label > 0, gray(256));
display_orientation(label, ori_map, 1, label > 0, gray(256));
for sigma = 2:10
    w_ori_smooth = imfilter(w_ori, fspecial('gaussian', 5*sigma, sigma));
    w_ori_map = 180*angle(w_ori_smooth)/pi;
    display_orientation(label, w_ori_map, 4, label > 0, gray(256));
end
%%
for sigma = 1:10
    display_weighted_orientation(line_map, ori_map, 1, [], sigma);
    saveas(gcf, ['C:\isbe\asymmetry_project\data\misc\figures\ori_map\map' zerostr(sigma,2) '.fig']);
    close(gcf);
    
    display_weighted_orientation(line_map, ori_map2, 1, [], sigma);
    saveas(gcf, ['C:\isbe\asymmetry_project\data\misc\figures\ori_map2\map' zerostr(sigma,2) '.fig']);
    close(gcf);
end
%%
for sigma = 1:10
%    
end
%%
display_weighted_orientation(line_map, ori_map, 1, line_map > 0, []); a1 = gca;
display_weighted_orientation(line_map, ori_map2, 1, line_map > 0, []); a2 = gca;
linkaxes([a1 a2]);
%%
display_weighted_orientation(line_map, ori_map, 1, [], []); a3 = gca;
display_weighted_orientation(line_map, ori_map2, 1, [], []); a4 = gca;
linkaxes([a3 a4]);
%%
for ii = 0:6
    sigma = 2^ii;
    ori_smooth = imfilter(w_ori, fspecial('gaussian', 5*sigma, sigma));
    abs_smooth = imfilter(abs(w_ori), fspecial('gaussian', 5*sigma, sigma));
    w_ori_smooth = abs_smooth .* exp(1i*angle(ori_smooth));
    figure;
    subplot(1,2,1); imagesc(abs_smooth); axis image;
    subplot(1,2,2); image(complex2rgb(w_ori_smooth.^2)); axis image;
end
%%
for ii = 0:6
    sigma = 2^ii;
    w_ori_smooth = imfilter(w_ori, fspecial('gaussian', 5*sigma, sigma));
    label = mb_non_maximal_supp(w_ori_smooth, ori_map, 1);
    display_weighted_orientation(w_ori_smooth, [], 1, label > 0, []);
end
%%
data_list = dir('C:\isbe\asymmetry_project\data\line_maps\2004_screening\normals\*CC*data.mat');
phi = linspace(0, pi, 100);

for ii = 1:2:39
    figure;
    subplot(3,2,5);
    line_map = imresize(u_load(['C:\isbe\asymmetry_project\data\line_maps\2004_screening\normals\' data_list(ii).name]), 0.5, 'bilinear');
    ori_map = imresize(u_load(['C:\isbe\asymmetry_project\data\orientation_maps\2004_screening\normals\' data_list(ii).name]), 0.5, 'bilinear');
    w_ori = line_map .* exp(1i * pi * ori_map / 180);
    [theta rho] = weighted_complex_rose(w_ori, 180);
    
    max_val = max(rho(:));
    
    %Plot radial axes marker
    plot([-max_val max_val], [0 0], 'k--'); hold on;
    plot([0 0], [0 max_val], 'k--');
    for jj = 2:2:10
        plot(jj*max_val*cos(phi)/10, jj*max_val*sin(phi)/10, 'k--');
    end

    plot(rho.*cos(theta), rho.*sin(theta)); axis equal;
    axis([-max_val max_val min(rho.*sin(theta)) max_val]);
    set(gca, 'xticklabel', []);
    title(['Orientation histogram for ' data_list(ii).name(1:3) ' ' data_list(ii).name(4:6)]);
        
    subplot(3,2,[1 3]); imagesc(line_map); axis image;
    clear line_map ori_map w_ori;
    
    subplot(3,2,6);
    line_map = imresize(u_load(['C:\isbe\asymmetry_project\data\line_maps\2004_screening\normals\' data_list(ii+1).name]), 0.5, 'bilinear');
    ori_map = imresize(u_load(['C:\isbe\asymmetry_project\data\orientation_maps\2004_screening\normals\' data_list(ii+1).name]), 0.5, 'bilinear');
    w_ori = line_map .* exp(1i * pi * ori_map / 180);
    [theta rho] = weighted_complex_rose(w_ori, 180);
    
    max_val = max(rho(:));
    
    %Plot radial axes marker
    plot([-max_val max_val], [0 0], 'k--'); hold on;
    plot([0 0], [0 max_val], 'k--');
    for jj = 2:2:10
        plot(jj*max_val*cos(phi)/10, jj*max_val*sin(phi)/10, 'k--');
    end

    plot(rho.*cos(theta), rho.*sin(theta)); axis equal;
    axis([-max_val max_val min(rho.*sin(theta)) max_val]);
    set(gca, 'xticklabel', []);
    title(['Orientation histogram for ' data_list(ii).name(1:3) ' ' data_list(ii).name(4:6)]);
    
    subplot(3,2,[2 4]); imagesc(line_map); axis image;
    clear line_map ori_map w_ori;
     
end
%%
mam_list = dir('C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\*.mat');
mam_names = get_mammo_info(mam_list);

mkdir C:\isbe\asymmetry_project\data\synthetic_backgrounds\film_bgs\
im = 1;
for ii = 1:length(mam_list);
    
    mammo = u_load(['C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\' mam_list(ii).name]);
    mask = u_load(['C:\isbe\asymmetry_project\data\masks\2004_screening\abnormals\' mam_list(ii).name(1:10) '_mask.mat']);
    
    if strcmpi(mam_names{ii}(4), 'R')
        mammo = fliplr(mammo);
        mask = fliplr(mask);
    end
    
    [r c] = size(mammo);
    r2 = round(r/2);
    
    c2 = find(any(mask), 1, 'last');
    
    if (c - c2) > 510
        bg = mammo(r2-255:r2+256, c2:c2+511);
        %figure; imagesc(bg); axis image;
        if strcmpi(mam_names{ii}(4), 'R')
            bg = fliplr(bg);
        end
        save(['C:\isbe\asymmetry_project\data\synthetic_backgrounds\film_bgs\bg' zerostr(im,3) '.mat'], 'bg');
        im = im+1;
    end
    clear mammo mask bg;
end
%%
rad_map = load_uint8('C:\isbe\asymmetry_project\data\weighted_radial_maps2\2004_screening\abnormals\002LCC_rad_map_016.mat');
tic;
for ii = 1:12
    rad_map(:,:,ii) = imfilter(rad_map(:,:,ii), fspecial('gaussian', 49, 8));
end
t2 = toc;
clear rad_map; pack;

rad_map = load_uint8('C:\isbe\asymmetry_project\data\weighted_radial_maps2\2004_screening\abnormals\002LCC_rad_map_016.mat');
tic;
rad_map = imfilter(rad_map, fspecial('gaussian', 25, 4));
t1 = toc;
clear rad_map; pack;
%%
%mkdir C:\isbe\asymmetry_project\data\ori_linop_maps\2004_screening\abnormals\
line_list = dir('C:\isbe\asymmetry_project\data\line_maps\2004_screening\abnormals\*.mat');
for ii = 1:length(line_list);
    display(['processing image ' num2str(ii)]);
    line_map = u_load(['C:\isbe\asymmetry_project\data\line_maps\2004_screening\abnormals\' line_list(ii).name]);
    [dummy ori_map] = line_operator_conv(line_map, 180, 11, 39, 'degrees');
    save(['C:\isbe\asymmetry_project\data\ori_linop_maps\2004_screening\abnormals\' line_list(ii).name], 'ori_map');
    
    clear ori_map dummy line_map;
end
    