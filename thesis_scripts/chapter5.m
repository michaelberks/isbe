%Chapter 5 script - produce results and figures for chapter 5 in thesis:
%separating a mass fromt eh background

%Examples we like are [16 26 28 39 44 57 58 59 64 75 104 106 108 121 133 177]
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 5.? in chapter
anno_list = dir('C:\isbe\dev\annotations\*.mat');

%%
cp1 = 1.0e+003*[0.8633 -2.1970 0.5504];
cp2 = 1.0e+003*[0.7854 -2.8339 0.6824];


for jj = 104%[16 26 28 39 44 57 58 59 64 75 104 106 108 121 133 177]
    anno = u_load(['C:\isbe\dev\annotations\', anno_list(jj).name]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Figure 5.? in chapter
    roi = cell(5,1);
    roi{1} = double(anno.mass_ROI);

    sigma = 18;
    gauss_filt = fspecial('gaussian', 5*sigma, sigma);
    roi{2} = imfilter(roi{1}, gauss_filt, 'symmetric');

    roi{3} = roi{2};
    roi{3}(anno.mass_list) = double(anno.mass_spline_it(anno.mass_list));

    roi{4} = double(anno.mass_spline_it);

    roi{5} = double(anno.mass_ROI - anno.mass_spline_it);

    c_lims = [min(roi{1}(:)) max(roi{1}(:))];
    z_lims = [min(roi{2}(:)) max(roi{2}(:))];
    %
    f1 = figure(...
        'Units', 'centimeters',...
        'Position', [0, 0, 14.4, 23.2],...
        'PaperPositionMode','auto',...
        'WindowStyle', 'normal',...
        'Color', [1, 1, 1],...
        'Visible', 'on');

    [r c] = size(roi{1}(50:650, 50:700));

    a_h = 4.5;
    a_w = a_h*c/r;
    buf = (7.2 - a_w)/3;

    for ii = 1:5
        axes_l(ii) = axes(...
        'Units', 'centimeters',...
        'position', [2*buf, (ii-1)*4.8, a_w, a_h]);

        imagesc(uint8(roi{6-ii}(50:650, 50:700))); colormap(gray(256)); 
        if ii>1, caxis(c_lims); end
        set(axes_l(ii), 'visible', 'off');

        axes_r(ii) = axes(...
        'Units', 'centimeters',...
        'position', [7.2, (ii-1)*4.8, 7.2, 4.8]);

        surf(roi{6-ii}(150:1:550, 150:1:600), 'EdgeColor', 'none'); colormap(gray(256));        
        if ii>1
            caxis(c_lims); set(axes_r(ii), 'Zlim', z_lims);
            set(axes_r(ii), 'units', 'normalized', 'CameraPosition', cp2);
        else
            set(axes_r(ii), 'units', 'normalized', 'CameraPosition', cp2);
        end
        set(axes_r(ii), 'visible', 'off');
    end
    % Now we need to put in some text
    
%     axes(axes_r(4)); text(0,1,'some writing here 2', 'units', 'normalized');
%     axes(axes_r(3)); text(0,1,'some writing here 3', 'units', 'normalized');
%     axes(axes_r(2)); text(0,1,'some writing here 4', 'units', 'normalized');
%     axes(axes_r(1)); text(0,1,'some writing here 5', 'units', 'normalized');

    %print('-dpdf', '-noui', '-painters', ['-f' num2str(f1)], '-r864', 'C:\isbe\thesis\figures\5\interp_new_method.pdf');
    %print('-deps2', '-noui', '-painters', ['-f' num2str(f1)], '-r300', 'C:\isbe\thesis\figures\5\interp_new_methoda.eps');
    print('-dtiff', '-noui', '-painters', ['-f' num2str(f1)], '-r300', 'C:\isbe\thesis\figures\5\interp_old_method.tif');
end
%%
cp1 = 1.0e+003*[0.8633 -2.1970 0.5504];
cp2 = 1.0e+003*[0.7854 -2.8339 0.6824];


for jj = 104%[16 26 28 39 44 57 58 59 64 75 104 106 108 121 133 177]
    anno = u_load(['C:\isbe\dev\annotations\', anno_list(jj).name]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Figure 5.? in chapter
    roi = cell(5,1);
    roi{1} = double(anno.mass_ROI);

    sigma = 18;
    gauss_filt = fspecial('gaussian', 5*sigma, sigma);
    roi{2} = imfilter(roi{1}, gauss_filt, 'symmetric');

    roi{3} = roi{2};
    roi{3}(anno.mass_list) = double(anno.mass_spline_it(anno.mass_list));

    roi{4} = double(anno.mass_sub_it);

    roi{5} = double(anno.mass_ROI) - anno.mass_sub_it;

    c_lims = [min(roi{1}(:)) max(roi{1}(:))];
    z_lims = [min(roi{2}(:)) max(roi{2}(:))];
    %
    f1 = figure(...
        'Units', 'centimeters',...
        'Position', [0, 0, 14.4, 23.2],...
        'PaperPositionMode','auto',...
        'WindowStyle', 'normal',...
        'Color', [1, 1, 1],...
        'Visible', 'on');

    [r c] = size(roi{1}(50:650, 50:700));

    a_h = 4.5;
    a_w = a_h*c/r;
    buf = (7.2 - a_w)/3;
    
    axes_r = zeros(5,1);
    axes_l = zeros(5,1);
    for ii = 1:5
        axes_l(ii) = axes(...
            'Units', 'centimeters',...
            'position', [2*buf, (ii-1)*4.8, a_w, a_h]);

        imagesc(uint8(roi{6-ii}(50:650, 50:700))); colormap(gray(256)); 
        if ~(ii==2), caxis(c_lims); end
        set(axes_l(ii), 'visible', 'off');
    end

    for ii = 1:5

        axes_r(ii) = axes(...
            'Units', 'centimeters',...
            'position', [7.2, (ii-1)*4.8, 7.2, 4.8]);

        surf(roi{6-ii}(150:1:550, 150:1:600), 'EdgeColor', 'none'); colormap(gray(256));        
        if ~(ii==2)
            caxis(c_lims); set(axes_r(ii), 'Zlim', z_lims);
            set(axes_r(ii), 'units', 'normalized', 'CameraPosition', cp2);
        else
            set(axes_r(ii), 'units', 'normalized', 'CameraPosition', cp2);
        end
        set(axes_r(ii), 'visible', 'off');
    end 
    
    %print('-dpdf', '-noui', '-painters', ['-f' num2str(f1)], '-r864', 'C:\isbe\thesis\figures\5\interp_new_method.pdf');
    %print('-deps2', '-noui', '-painters', ['-f' num2str(f1)], '-r300', 'C:\isbe\thesis\figures\5\interp_new_methoda.eps');
    print('-dtiff', '-noui', '-painters', ['-f' num2str(f1)], '-r300', 'C:\isbe\thesis\figures\5\interp_new_method.tif');
end
%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
anno1 = u_load(['C:\isbe\dev\annotations\', anno_list(16).name]);
anno2 = u_load(['C:\isbe\dev\annotations\', anno_list(26).name]);
anno3 = u_load(['C:\isbe\dev\annotations\', anno_list(104).name]);
anno4 = u_load(['C:\isbe\dev\annotations\', anno_list(106).name]);

f1 = figure(...
        'Units', 'centimeters',...
        'Position', [0, 0, 14.4, 19.4],...
        'PaperPositionMode','auto',...
        'WindowStyle', 'normal',...
        'Color', [1, 1, 1],...
        'Visible', 'on');
colormap(gray(256));

axes('Units', 'centimeters', 'position', [0, 0, 4.7, 4.7]);
imagesc(anno1.mass_ROI); set(gca, 'visible', 'off');
axes('Units', 'centimeters', 'position', [4.85, 0, 4.7, 4.7]); 
imagesc(uint8(double(anno1.mass_ROI)-anno1.mass_sub_it)); set(gca, 'visible', 'off');
axes('Units', 'centimeters', 'position', [9.7, 0, 4.7, 4.7]);
imagesc(uint8(anno1.mass_sub_it)); set(gca, 'visible', 'off');

axes('Units', 'centimeters', 'position', [0, 4.9, 4.7, 4.7]);
imagesc(anno2.mass_ROI); set(gca, 'visible', 'off');
axes('Units', 'centimeters', 'position', [4.85, 4.9, 4.7, 4.7]);
imagesc(uint8(double(anno2.mass_ROI)-anno2.mass_sub_it)); set(gca, 'visible', 'off');
axes('Units', 'centimeters', 'position', [9.7, 4.9, 4.7, 4.7]);
imagesc(uint8(anno2.mass_sub_it)); set(gca, 'visible', 'off');

axes('Units', 'centimeters', 'position', [0, 9.8, 4.7, 4.7]);
imagesc(anno3.mass_ROI); set(gca, 'visible', 'off');
axes('Units', 'centimeters', 'position', [4.85, 9.8, 4.7, 4.7]);
imagesc(uint8(double(anno3.mass_ROI)-anno3.mass_sub_it)); set(gca, 'visible', 'off');
axes('Units', 'centimeters', 'position', [9.7, 9.8, 4.7, 4.7]);
imagesc(uint8(anno3.mass_sub_it)); set(gca, 'visible', 'off');

axes('Units', 'centimeters', 'position', [0, 14.7, 4.7, 4.7]);
imagesc(anno4.mass_ROI); set(gca, 'visible', 'off');
axes('Units', 'centimeters', 'position', [4.85, 14.7, 4.7, 4.7]);
imagesc(uint8(double(anno4.mass_ROI)-anno4.mass_sub_it)); set(gca, 'visible', 'off');
axes('Units', 'centimeters', 'position', [9.7, 14.7, 4.7, 4.7]);
imagesc(uint8(anno4.mass_sub_it)); set(gca, 'visible', 'off');

print('-dtiff', '-noui', '-painters', ['-f' num2str(f1)], '-r300', 'C:\isbe\thesis\figures\5\interpolation_new_masses_16_26_104_106.tif'); 
%%
anno1 = u_load(['C:\isbe\dev\annotations\', anno_list(16).name]);
anno2 = u_load(['C:\isbe\dev\annotations\', anno_list(26).name]);
anno3 = u_load(['C:\isbe\dev\annotations\', anno_list(104).name]);
anno4 = u_load(['C:\isbe\dev\annotations\', anno_list(106).name]);

f1 = figure(...
        'Units', 'centimeters',...
        'Position', [0, 0, 14.4, 19.4],...
        'PaperPositionMode','auto',...
        'WindowStyle', 'normal',...
        'Color', [1, 1, 1],...
        'Visible', 'on');
colormap(gray(256));

axes('Units', 'centimeters', 'position', [0, 0, 4.7, 4.7]);
imagesc(anno1.mass_ROI); set(gca, 'visible', 'off');
axes('Units', 'centimeters', 'position', [4.85, 0, 4.7, 4.7]); 
imagesc(uint8(anno1.mass_spline_it)); set(gca, 'visible', 'off');
axes('Units', 'centimeters', 'position', [9.7, 0, 4.7, 4.7]);
imagesc(uint8(anno1.mass_sub)); set(gca, 'visible', 'off');

axes('Units', 'centimeters', 'position', [0, 4.9, 4.7, 4.7]);
imagesc(anno2.mass_ROI); set(gca, 'visible', 'off');
axes('Units', 'centimeters', 'position', [4.85, 4.9, 4.7, 4.7]);
imagesc(uint8(anno2.mass_spline_it)); set(gca, 'visible', 'off');
axes('Units', 'centimeters', 'position', [9.7, 4.9, 4.7, 4.7]);
imagesc(uint8(anno2.mass_sub)); set(gca, 'visible', 'off');

axes('Units', 'centimeters', 'position', [0, 9.8, 4.7, 4.7]);
imagesc(anno3.mass_ROI); set(gca, 'visible', 'off');
axes('Units', 'centimeters', 'position', [4.85, 9.8, 4.7, 4.7]);
imagesc(uint8(anno3.mass_spline_it)); set(gca, 'visible', 'off');
axes('Units', 'centimeters', 'position', [9.7, 9.8, 4.7, 4.7]);
imagesc(uint8(anno3.mass_sub)); set(gca, 'visible', 'off');

axes('Units', 'centimeters', 'position', [0, 14.7, 4.7, 4.7]);
imagesc(anno4.mass_ROI); set(gca, 'visible', 'off');
axes('Units', 'centimeters', 'position', [4.85, 14.7, 4.7, 4.7]);
imagesc(uint8(anno4.mass_spline_it)); set(gca, 'visible', 'off');
axes('Units', 'centimeters', 'position', [9.7, 14.7, 4.7, 4.7]);
imagesc(uint8(anno4.mass_sub)); set(gca, 'visible', 'off');

print('-dtiff', '-noui', '-painters', ['-f' num2str(f1)], '-r300', 'C:\isbe\thesis\figures\5\interpolation_old_masses_16_26_104_106.tif'); 
%%
%--------------------------------------------------------------------------
anno1 = u_load(['C:\isbe\dev\annotations\', anno_list(16).name]);
anno2 = u_load(['C:\isbe\dev\annotations\', anno_list(26).name]);
anno3 = u_load(['C:\isbe\dev\annotations\', anno_list(104).name]);
anno4 = u_load(['C:\isbe\dev\annotations\', anno_list(106).name]);

figure; imagesc(anno1.mass_ROI); axis image; colormap(gray(256)); hold on;
plot(anno1.mass_outline(:,1), anno1.mass_outline(:,2), 'b');
for ii = 1:length(anno1.mass_spicules)
    plot(anno1.mass_spicules(ii).outline(:,1), anno1.mass_spicules(ii).outline(:,2), 'r');
end
figure; imagesc(anno4.mass_ROI); axis image; colormap(gray(256)); hold on;
plot(anno4.mass_outline(:,1), anno4.mass_outline(:,2), 'b');
for ii = 1:length(anno4.mass_spicules)
    plot(anno4.mass_spicules(ii).outline(:,1), anno4.mass_spicules(ii).outline(:,2), 'r');
end
figure; imagesc(anno2.mass_ROI); axis image; colormap(gray(256)); hold on;
plot(anno2.mass_outline(:,1), anno2.mass_outline(:,2), 'b');
for ii = 1:length(anno2.mass_spicules)
    plot(anno2.mass_spicules(ii).outline(:,1), anno2.mass_spicules(ii).outline(:,2), 'r');
end
figure; imagesc(anno3.mass_ROI); axis image; colormap(gray(256)); hold on;
plot(anno3.mass_outline(:,1), anno3.mass_outline(:,2), 'b');
for ii = 1:length(anno3.mass_spicules)
    plot(anno3.mass_spicules(ii).outline(:,1), anno3.mass_spicules(ii).outline(:,2), 'r');
end
%%
roi = anno1.mass_ROI(98:end-126,:);
[r c] = size(roi);
f1 = figure(...
        'Units', 'pixels',...
        'Position', [0, 0, r/2, c/2],...
        'PaperPositionMode','auto',...
        'WindowStyle', 'normal',...
        'Color', [1, 1, 1],...
        'Visible', 'on');

axes('Units', 'pixels', 'position', [0, 0, r/2, c/2]);
imagesc(roi); colormap(gray(256)); hold on;
plot(anno1.mass_outline(:,1), anno1.mass_outline(:,2)-97, 'b', 'LineWidth', 3.0);
for ii = 1:length(anno4.mass_spicules)
    plot(anno1.mass_spicules(ii).outline(:,1), anno1.mass_spicules(ii).outline(:,2)-97, 'r', 'LineWidth', 3.0);
end
set(gca, 'visible', 'off');
print('-dtiff', '-noui', '-painters', ['-f' num2str(f1)], '-r300', 'C:\isbe\thesis\figures\5\mass_16_annotated.tif');

roi = anno2.mass_ROI(1:end-20,:);
[r c] = size(roi);
f1 = figure(...
        'Units', 'pixels',...
        'Position', [0, 0, r/2, c/2],...
        'PaperPositionMode','auto',...
        'WindowStyle', 'normal',...
        'Color', [1, 1, 1],...
        'Visible', 'on');

axes('Units', 'pixels', 'position', [0, 0, r/2, c/2]);
imagesc(roi); colormap(gray(256)); hold on;
plot(anno2.mass_outline(:,1), anno2.mass_outline(:,2), 'b', 'LineWidth', 3.0);
for ii = 1:length(anno2.mass_spicules)
    plot(anno2.mass_spicules(ii).outline(:,1), anno2.mass_spicules(ii).outline(:,2), 'r', 'LineWidth', 3.0);
end
set(gca, 'visible', 'off');
print('-dtiff', '-noui', '-painters', ['-f' num2str(f1)], '-r300', 'C:\isbe\thesis\figures\5\mass_26_annotated.tif');

roi = anno3.mass_ROI(44:end,:);
[r c] = size(roi);
f1 = figure(...
        'Units', 'pixels',...
        'Position', [0, 0, r/2, c/2],...
        'PaperPositionMode','auto',...
        'WindowStyle', 'normal',...
        'Color', [1, 1, 1],...
        'Visible', 'on');

axes('Units', 'pixels', 'position', [0, 0, r/2, c/2]);
imagesc(roi); colormap(gray(256)); hold on;
plot(anno3.mass_outline(:,1), anno3.mass_outline(:,2)-43, 'b', 'LineWidth', 3.0);
for ii = 1:length(anno3.mass_spicules)
    plot(anno3.mass_spicules(ii).outline(:,1), anno3.mass_spicules(ii).outline(:,2)-43, 'r', 'LineWidth', 3.0);
end
set(gca, 'visible', 'off');
print('-dtiff', '-noui', '-painters', ['-f' num2str(f1)], '-r300', 'C:\isbe\thesis\figures\5\mass_104_annotated.tif');

roi = anno4.mass_ROI(:,136:end-58);
[r c] = size(roi);
f1 = figure(...
        'Units', 'pixels',...
        'Position', [0, 0, r/2, c/2],...
        'PaperPositionMode','auto',...
        'WindowStyle', 'normal',...
        'Color', [1, 1, 1],...
        'Visible', 'on');

axes('Units', 'pixels', 'position', [0, 0, r/2, c/2]);
imagesc(roi); colormap(gray(256)); hold on;
plot(anno4.mass_outline(:,1)-135, anno4.mass_outline(:,2), 'b', 'LineWidth', 3.0);
for ii = 1:length(anno4.mass_spicules)
    plot(anno4.mass_spicules(ii).outline(:,1)-135, anno4.mass_spicules(ii).outline(:,2), 'r', 'LineWidth', 3.0);
end
set(gca, 'visible', 'off');
print('-dtiff', '-noui', '-painters', ['-f' num2str(f1)], '-r300', 'C:\isbe\thesis\figures\5\mass_106_annotated.tif');
%%
write_im_from_colormap(anno1.mass_ROI(98:end-126,:), 'C:\isbe\thesis\figures\5\mass_roi_16.bmp', colormap(gray(256)));
write_im_from_colormap(uint8(anno1.mass_sub_it(98:end-126,:)), 'C:\isbe\thesis\figures\5\mass_sub_16.bmp', colormap(gray(256)));
write_im_from_colormap(double(anno1.mass_ROI(98:end-126,:))-anno1.mass_sub_it(98:end-126,:), 'C:\isbe\thesis\figures\5\mass_bg_16.bmp', colormap(gray(256)));

write_im_from_colormap(anno2.mass_ROI(1:end-20,:), 'C:\isbe\thesis\figures\5\mass_roi_26.bmp', colormap(gray(256)));
write_im_from_colormap(uint8(anno2.mass_sub_it(1:end-20,:)), 'C:\isbe\thesis\figures\5\mass_sub_26.bmp', colormap(gray(256)));
write_im_from_colormap(double(anno2.mass_ROI(1:end-20,:))-anno2.mass_sub_it(1:end-20,:), 'C:\isbe\thesis\figures\5\mass_bg_26.bmp', colormap(gray(256)));

write_im_from_colormap(anno3.mass_ROI(44:end,:), 'C:\isbe\thesis\figures\5\mass_roi_104.bmp', colormap(gray(256)));
write_im_from_colormap(uint8(anno3.mass_sub_it(44:end,:)), 'C:\isbe\thesis\figures\5\mass_sub_104.bmp', colormap(gray(256)));
write_im_from_colormap(double(anno3.mass_ROI(44:end,:))-anno3.mass_sub_it(44:end,:), 'C:\isbe\thesis\figures\5\mass_bg_104.bmp', colormap(gray(256)));

write_im_from_colormap(anno4.mass_ROI(:,136:end-58), 'C:\isbe\thesis\figures\5\mass_roi_106.bmp', colormap(gray(256)));
write_im_from_colormap(uint8(anno4.mass_sub_it(:,136:end-58)), 'C:\isbe\thesis\figures\5\mass_sub_106.bmp', colormap(gray(256)));
write_im_from_colormap(double(anno4.mass_ROI(:,136:end-58))-anno4.mass_sub_it(:,136:end-58), 'C:\isbe\thesis\figures\5\mass_bg_106.bmp', colormap(gray(256)));


%--------------------------------------------------------------------------
%**************************************************************************
%%
%First up, look at all the mass regions, and the final separations to pick
%some good example figure
mass_list = dir('C:\isbe\dev\masses1024x1024\*.mat');
%%
for ii = 151:180
    load(['C:\isbe\dev\masses1024x1024\', mass_list(ii).name]);
    figure;
    subplot(1,3,1); imagesc(mass.mass_ROI); axis image; colormap(gray(256));
    subplot(1,3,2); imagesc(mass.background_ROI); axis image; colormap(gray(256));
    subplot(1,3,3); imagesc(uint8(mass.mass_ROI-mass.background_ROI)); axis image; colormap(gray(256));
    clear mass;
end
%%
for ii = [16 26 28 39 44 57 58 59 64 75 104 106 108 121 133 177]
    load(['C:\isbe\dev\annotations\', anno_list(ii).name]);
    figure; surf(double(mass.mass_ROI(1:4:end, 1:4:end)), 'EdgeColor', 'none'); axis image; colormap(gray(256));
    clear mass;
end
%%
%Examples we like are
for ii = [16 26 28 39 44 57 58 59 64 75 104 106 108 121 133 177]
    load(['C:\isbe\dev\annotations\', anno_list(ii).name]);
    figure; imagesc(mass.mass_ROI); axis image; colormap(gray(256));
    clear mass;
end
%%
for ii = [16 26 28 39 44 57 58 59 64 75 104 106 108 121 133 177]
    load(['C:\isbe\dev\annotations\', anno_list(ii).name]);
    figure;
    subplot(1,3,1); imagesc(mass.mass_ROI); axis image; colormap(gray(256));
    subplot(1,3,2); imagesc(uint8(mass.mass_sub_it)); axis image; colormap(gray(256));
    subplot(1,3,3); imagesc(uint8(double(mass.mass_ROI)-mass.mass_sub_it)); axis image; colormap(gray(256));
    clear mass;
end
%%
%39 could be a winner baby - that is m04_034LMLa. WE now want to load the
%annotation object for this mass
anno = u_load('C:\isbe\dev\annotations\an04_034LMLa.mat');
%%
%Mass 141 is good for showing the profile example
anno_list = dir('C:\isbe\dev\annotations\*.mat');
anno = u_load(['C:\isbe\dev\annotations\', anno_list(140).name]);
for ii = 347
    figure; hold on;
    plot(anno.mass_ROI(ii,:)');
    plot(anno.mass_spline(ii,:)', 'r:');
    plot(anno.mass_spline_it(ii,:)', 'g:');
end
%
n1 = 50; n2 = 10; spacing = 10; sigma = 18; plot_flag = 1; green = 'biharmTPS'; f_method = 1;
gauss_filt = fspecial('gaussian', 5*sigma, sigma);
shape_bw = roipoly(anno.mass_ROI, anno.mass_outline(:,1), anno.mass_outline(:,2));
%
figure; hold on;
mass_ROI = anno.mass_ROI;
smooth_ROI_0 = imfilter(mass_ROI, gauss_filt, 'symmetric');
plot(mass_ROI(347,:)');

[bg_estimates_0 p_list] = spline_estimation(smooth_ROI_0, shape_bw, n1, n2, spacing, green);
mass_ROI(p_list) = uint8(bg_estimates_0);
smooth_ROI_1 = imfilter(mass_ROI, gauss_filt, 'symmetric');
plot(mass_ROI(347,:)', 'r:');

[bg_estimates_1] = spline_estimation(smooth_ROI_1, shape_bw, n1, n2, spacing, green);
mass_ROI(p_list) = uint8(bg_estimates_1);
smooth_ROI_2 = imfilter(mass_ROI, gauss_filt, 'symmetric');
plot(mass_ROI(347,:)', 'g:');

[bg_estimates_2] = spline_estimation(smooth_ROI_2, shape_bw, n1, n2, spacing, green);
mass_ROI(p_list) = uint8(bg_estimates_2);
smooth_ROI_3 = imfilter(mass_ROI, gauss_filt, 'symmetric');
plot(mass_ROI(347,:)', 'm:');

[bg_estimates_3] = spline_estimation(smooth_ROI_3, shape_bw, n1, n2, spacing, green);
mass_ROI(p_list) = uint8(bg_estimates_3);
smooth_ROI_4 = imfilter(mass_ROI, gauss_filt, 'symmetric');
plot(mass_ROI(347,:)', 'c:');
%
smooth_ROI_0a = smooth_ROI_0;
smooth_ROI_0a(p_list) = bg_estimates_0;

smooth_ROI_1a = smooth_ROI_1;
smooth_ROI_1a(p_list) = bg_estimates_1;

smooth_ROI_2a = smooth_ROI_3;
smooth_ROI_2a(p_list) = bg_estimates_2;

smooth_ROI_3a = smooth_ROI_3;
smooth_ROI_3a(p_list) = bg_estimates_3;
%%
f1 = figure(...
    'Units', 'centimeters',...
    'Position', [0, 0, 10, 8],...
    'PaperPositionMode','auto',...
    'WindowStyle', 'normal',...
    'Color', [1, 1, 1],...
    'Visible', 'on'); 
hold on;
plot(smooth_ROI_0(347,:)', 'b:', 'LineWidth', 1.5);
plot(smooth_ROI_0a(347,:)', 'r:', 'LineWidth', 1.5);
plot(smooth_ROI_1a(347,:)', 'm:', 'LineWidth', 1.5);
plot(smooth_ROI_2a(347,:)', 'g:', 'LineWidth', 1.5);
plot(smooth_ROI_3a(347,:)', 'b', 'LineWidth', 1.5);
ylabel('Grey-level intensity');
xlabel('Distance along profile (in pixels)');
%print_eps('C:\isbe\thesis\figures\5\iteration_profiles.eps');
print('-dtiff', '-noui', '-painters', f1, '-r300', 'C:\isbe\thesis\figures\5\iteration_profiles.tif');
%
f1 = figure(...
    'Units', 'centimeters',...
    'Position', [0, 0, 10, 8],...
    'PaperPositionMode','auto',...
    'WindowStyle', 'normal',...
    'Color', [1, 1, 1],...
    'Visible', 'on'); 
hold on;
plot(smooth_ROI_0a(347,:)', 'r:', 'LineWidth', 1.5);
plot(smooth_ROI_0(347,:)', 'b', 'LineWidth', 1.5);
ylabel('Grey-level intensity');
xlabel('Distance along profile (in pixels)');
%print_eps('C:\isbe\thesis\figures\5\estimation_profiles.eps');
print('-dtiff', '-noui', '-painters', f1, '-r300', 'C:\isbe\thesis\figures\5\estimation_profiles.tif');
%%
[r c] = size(anno.mass_ROI);
f1 = figure(...
        'Units', 'pixels',...
        'Position', [0, 0, r/4, c/4],...
        'PaperPositionMode','auto',...
        'WindowStyle', 'normal',...
        'Color', [1, 1, 1],...
        'Visible', 'on');

axes('Units', 'pixels', 'position', [0, 0, r/4, c/4]);
imagesc(anno.mass_ROI); colormap(gray(256)); hold on;
plot([0.5 size(anno.mass_ROI,2)], [347 347], 'y:', 'LineWidth', 2.0);
set(gca, 'visible', 'off');
print('-dtiff', '-noui', '-painters', ['-f' num2str(f1)], '-r300', 'C:\isbe\thesis\figures\5\mass_roi_profile.tif');

%%
for ii = 200:50:600
    figure; hold on;
    plot(smooth_ROI_0(:,ii));
    plot(smooth_ROI_0a(:,ii), 'b:');
    plot(smooth_ROI_1a(:,ii), 'r:');
    plot(smooth_ROI_2a(:,ii), 'g:');
    plot(smooth_ROI_3a(:,ii), 'm:');
end
%%

figure; imagesc(smooth_ROI_0a); axis image; colormap(gray(256));
figure; imagesc(smooth_ROI_3a); axis image; colormap(gray(256));

figure; surf(double(smooth_ROI_0a(1:4:end, 1:4:end))); colormap(gray(256));
figure; surf(double(smooth_ROI_3a(1:4:end, 1:4:end))); colormap(gray(256));
%%
figure; imagesc(smooth_ROI_0); axis image; colormap(gray(256));
figure; imagesc(smooth_ROI_1); axis image; colormap(gray(256));
figure; imagesc(smooth_ROI_2); axis image; colormap(gray(256));
figure; imagesc(smooth_ROI_3); axis image; colormap(gray(256));
figure; imagesc(smooth_ROI_4); axis image; colormap(gray(256));
%%
figure; hold on;
plot(smooth_ROI_0(ii,:)');
plot(smooth_ROI_1(ii,:)', 'r:');
plot(smooth_ROI_2(ii,:)', 'g:');
plot(smooth_ROI_3(ii,:)', 'm:');
plot(smooth_ROI_4(ii,:)', 'c:');
%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%Make test data
mass_list = dir('C:\isbe\dev\annotations\*.mat');
normal_list = dir('C:\isbe\dev\background\images\normal1024\*.bmp');
r_mass = randperm(length(mass_list));
r_normal = randperm(length(normal_list));

%
mkdir C:\isbe\thesis\figures\5\experiment

mm = 1; nn = 1;
while nn <= 10
    bg = imread(['C:\isbe\dev\background\images\normal1024\', normal_list(r_normal(nn)).name]);
    mass = u_load(['C:\isbe\dev\annotations\', mass_list(r_mass(mm)).name]);
    [r c] = size(mass.mass_ROI);
    if r <= 1024 && c < 1024
        bg = double(bg(1:r, 1:c));
        test_data.mass_ROI = bg + mass.mass_sub_it;
        test_data.mass_outline = mass.mass_outline;
        mass_bw = roipoly(bg, mass.mass_outline(:,1), mass.mass_outline(:,2));
        test_data.real_bg = bg(mass_bw);
        save(['C:\isbe\thesis\figures\5\experiment\test_mass', zerostr(nn,2)], 'test_data');
        clear test_data
        nn = nn+1;
    end
    
    mm = mm+1;
end


%%
test_list = dir('C:\isbe\thesis\figures\5\experiment\*.mat');

[n1 sigma] = meshgrid(0:10:50, 5:5:30);
errors_original = zeros(size(n1));
errors_iterative = zeros(size(n1));
for ii = 1:numel(n1)
    [errors] = subtract_mass_test(test_list,'C:\isbe\thesis\figures\5\experiment\', n1(ii), 20, 10, sigma(ii), 1, 'biharmTPS');

    errors_original(ii) = mean(errors(:,1));
    errors_iterative(ii) = mean(errors(:,2));
end

%original best
[min_original i_o] = min(errors_original(:));
display(['Best original parameters are n1 = ', num2str(n1(i_o)), ' sigma = ', num2str(sigma(i_o))]);

%best iterative is:
[min_iterative i_i] = min(errors_iterative(:));
display(['Best iterative parameters are n1 = ', num2str(n1(i_i)), ' sigma = ', num2str(sigma(i_i))]);
%%
load('C:\isbe\thesis\figures\5\experiment\errors.mat');

f1 = figure(...
        'Units', 'pixels',...
        'Position', [10, 10, 600, 450],...
        'PaperPositionMode','auto',...
        'WindowStyle', 'normal',...
        'Color', [1, 1, 1],...
        'Visible', 'on');
plot(errors_original', 'LineWidth', 1.5);
set(gca, 'FontSize', 14, 'xlim', [1 11], 'xtick', 1:11, 'xticklabel', (0:10:100));
xlabel('n_1');
ylabel('Mean background estimation error');
legend({'\sigma = 5', '\sigma = 10', '\sigma = 15', '\sigma = 25', '\sigma = 25', '\sigma = 30', '\sigma = 35'});
print('-dtiff', '-noui', '-painters', gcf, '-r300', 'C:\isbe\thesis\figures\5\parameter_experiment.tif');
%%
mass_list = dir('C:\isbe\dev\masses1024x1024\*.mat');
for ii = 31:60
    load(['C:\isbe\dev\masses1024x1024\', mass_list(ii).name]);
    figure; imagesc(mass.background_ROI); colormap(gray(256)); axis image;
end
    
    