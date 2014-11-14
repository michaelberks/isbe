% Intsruction used to generate figures for ECR poster

%Get images for modifying a region and inserting a new mass
load('C:\isbe\dev\background\misc\mammogram_region_modified_mass_added.mat')

mammo = double(imread('C:\isbe\mammograms\new_CAD\BMP_2004\o04_008LCC.bmp'));
mammo_patch = mammo(mass.R1:mass.R2, mass.C1:mass.C2); clear mammo;
mammo_patch_mass = mammo_patch+mass.subtract_ROI;
modified_patch = mass.background_ROI;
modified_patch_mass = mass.background_ROI+mass.subtract_ROI;
new_mass = mass.subtract_ROI;

figure; imagesc(new_mass); axis image; colormap(gray(256));
figure; imagesc(mammo_patch); axis image; colormap(gray(256)); caxis([100 180]);
figure; imagesc(modified_patch); axis image; colormap(gray(256)); caxis([100 180]);
figure; imagesc(mammo_patch_mass); axis image; colormap(gray(256)); caxis([100 180]);
figure; imagesc(modified_patch_mass); axis image; colormap(gray(256)); caxis([100 180]);

write_im_from_colormap(new_mass, 'K:\isbe\conferences_and_symposia\ecr\new_mass.bmp', gray(256));
write_im_from_colormap(mammo_patch, 'K:\isbe\conferences_and_symposia\ecr\mammo_patch.bmp', gray(256), [100 180]);
write_im_from_colormap(mammo_patch_mass, 'K:\isbe\conferences_and_symposia\ecr\mammo_patch_mass.bmp', gray(256), [100 180]);
write_im_from_colormap(modified_patch, 'K:\isbe\conferences_and_symposia\ecr\modified_patch.bmp', gray(256), [100 180]);
write_im_from_colormap(modified_patch_mass, 'K:\isbe\conferences_and_symposia\ecr\modified_patch_mass.bmp', gray(256), [100 180]);

% Make animated GIF of adding a mass to a un/modified region 
figure('Position', [100,100, 662, 348], 'WindowStyle', 'normal', 'Color', [1, 1, 1]);
a1 = axes('Units', 'pixels', 'position', [20, 20, 301, 308]);
a2 = axes('Units', 'pixels', 'position', [341, 20, 301, 308]);
axes(a1); imagesc(mammo_patch); axis image; colormap(gray(256)); caxis([100 180]); axis off;
axes(a2); imagesc(modified_patch); axis image; colormap(gray(256)); caxis([100 180]); axis off;
frame1 = getframe(gcf);
gif1 = frame2im(frame1);
[gif1a map] = rgb2ind(gif1, 2^16);
imwrite(gif1a, map, 'K:\isbe\conferences_and_symposia\ecr\figures\adding_a_mass.gif', 'gif', 'WriteMode', 'overwrite', 'Loop', Inf, 'DelayTime', 1.0);
figure('Position', [100,100, 662, 348], 'WindowStyle', 'normal', 'Color', [1, 1, 1]);
a1 = axes('Units', 'pixels', 'position', [20, 20, 301, 308]);
a2 = axes('Units', 'pixels', 'position', [341, 20, 301, 308]);
axes(a1); imagesc(mammo_patch_mass); axis image; colormap(gray(256)); caxis([100 180]); axis off;
axes(a2); imagesc(modified_patch_mass); axis image; colormap(gray(256)); caxis([100 180]); axis off;
frame2 = getframe(gcf);
gif2 = frame2im(frame2);
[gif2a map] = rgb2ind(gif2, 2^16);
imwrite(gif2a, map, 'K:\isbe\conferences_and_symposia\ecr\figures\adding_a_mass.gif', 'gif', 'WriteMode', 'append', 'DelayTime', 1.0);
%%
clear
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make animated GIFs for real or synthetic masses
for m = 1:6
    mass_name = ['K:\isbe\conferences_and_symposia\iwdm2008\pyramid\poster\figures\poster_masses\syn_mass1_', zerostr(m,2), '.bmp'];
    gif_name = ['K:\isbe\conferences_and_symposia\ecr\figures\real_or_synthetic_s', zerostr(m,2), '.gif'];
    syn_mass = imread(mass_name);

    f1 = figure('Position', [100,100, 756, 756], 'WindowStyle', 'normal', 'Color', [1, 1, 1]);
    a1 = axes('Units', 'pixels', 'position', [0, 0, 756, 756]);
    imagesc(rgb2gray(syn_mass)); axis image; axis off; colormap(gray(256));

    frame1 = getframe(f1);
    gif1 = frame2im(frame1);
    [gif1a map] = rgb2ind(gif1, 2^16);
    imwrite(gif1a, map, gif_name, 'gif', 'WriteMode', 'overwrite', 'Loop', Inf, 'DelayTime', 5.0);
    
    axes(a1); text(20, 40, 'SYNTHETIC', 'FontSize', 36, 'Color', 'r');
    frame2 = getframe(f1);
    gif2 = frame2im(frame2);
    [gif2a map] = rgb2ind(gif2, 2^16);
    imwrite(gif2a, map, gif_name, 'WriteMode', 'append', 'DelayTime', 2.0);
    
    close(f1);
end
%%
for m = 1:6
    mass_name = ['K:\isbe\conferences_and_symposia\iwdm2008\pyramid\poster\figures\poster_masses\real_mass', zerostr(m,2), '.bmp'];
    gif_name = ['K:\isbe\conferences_and_symposia\ecr\figures\real_or_synthetic_r', zerostr(m,2), '.gif'];
    syn_mass = imread(mass_name);

    f1 = figure('Position', [100,100, 756, 756], 'WindowStyle', 'normal', 'Color', [1, 1, 1]);
    a1 = axes('Units', 'pixels', 'position', [0, 0, 756, 756]);
    imagesc(rgb2gray(syn_mass)); axis image; axis off; colormap(gray(256));

    frame1 = getframe(f1);
    gif1 = frame2im(frame1);
    [gif1a map] = rgb2ind(gif1, 2^16);
    imwrite(gif1a, map, gif_name, 'gif', 'WriteMode', 'overwrite', 'Loop', Inf, 'DelayTime', 5.0);
    
    axes(a1); text(20, 40, 'REAL', 'FontSize', 36, 'Color', 'g');
    frame2 = getframe(f1);
    gif2 = frame2im(frame2);
    [gif2a map] = rgb2ind(gif2, 2^16);
    imwrite(gif2a, map, gif_name, 'WriteMode', 'append', 'DelayTime', 2.0);
    
    close(f1);
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
test_im = randn(256);

f1 = figure('Position', [100,100, 400, 400], 'WindowStyle', 'normal', 'Color', [0, 0, 0]); 
imagesc(test_im); axis image; colormap(gray(256)); hold on;
plot(0:256, 0:256, 'b', 'LineWidth', 2);
test_frames(1) = getframe(f1);
plot(0:8:256, 0:8:256, 'rx', 'MarkerSize', 10);
test_frames(2) = getframe(f1);
close(f1);

for f = 1:length(test_frames)
    if f == 1
        write_mode = 'overwrite';
    else
        write_mode = 'append';
    end
    [X] = frame2im(test_frames(f));
    [gif_im, cmap] = rgb2ind(X, 65536);
    imwrite(gif_im, cmap, 'C:\isbe\test_movie.gif', 'gif', 'WriteMode', write_mode);
end
%test_frames(2) = getframe(gcf);

%movie2avi(test_frames, 'C:\isbe\test_movie.avi','compression','none','fps',10);