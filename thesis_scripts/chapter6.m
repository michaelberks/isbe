%Chapter 6 script
anno_list = dir('C:\isbe\dev\annotations\*.mat');
mass_list = dir('C:\isbe\dev\masses\*.mat');
u_list = u_load('C:\isbe\dev\files\u_files.mat');
%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Figure showing how points are allocated on a mass border

n_pts = 50;

mass = u_load(['C:\isbe\dev\masses\', u_list(7).name]);
mass_outline = mass.mass_outline;
mass_centroid = mass.mass_centroid;
nipple = mass.nipple;

c_to_n = (nipple - mass_centroid) / sqrt(sum((nipple - mass_centroid).^2));
c_to_n = repmat(c_to_n, size(mass_outline, 1), 1);

c_to_o = mass_outline - repmat(mass_centroid, size(mass_outline,1), 1);
c_to_o = c_to_o ./ repmat(sqrt(sum(c_to_o'.^2))', 1, 2);

[dummy m_idx1] = min(sum((c_to_o - c_to_n)'.^2));
[dummy m_idx2] = min(sum((c_to_o + c_to_n)'.^2));

dist = sqrt(sum((mass_centroid - mass_outline(1,:)).^2));
p1 = mass.mass_centroid + 1.5*dist*c_to_n(1,:);
p2 = mass.mass_centroid - dist*c_to_n(1,:);

idx1 = round(linspace(1, m_idx2, n_pts/2 + 1));
idx2 = round(linspace(m_idx2, size(mass_outline,1), n_pts/2 + 1));



x_min = min(mass_outline(:,1))-100;
x_max = max(mass_outline(:,1))+200;

y_min = min(mass_outline(:,2))-50;
y_max = max(mass_outline(:,2))+50;

f1 = figure(...
    'Units', 'pixels',...
    'Position', [0 0 x_max-x_min y_max-y_min],...
    'PaperPositionMode','auto',...
    'WindowStyle', 'normal',...
    'Color', [1, 1, 1],...
    'Visible', 'on');

a1 = axes(...
    'Units', 'pixels',...
    'Position', [0 0 x_max-x_min y_max-y_min],...
    'XLim', [x_min x_max],...
    'YLim', [y_min y_max],...
    'NextPlot', 'add',...
    'Visible', 'off');
    
plot(mass_outline(:,1), mass_outline(:,2), 'k', 'LineWidth', 3);
plot(mass_outline(idx1(1:end-1),1), mass_outline(idx1(1:end-1),2), 'r.', 'MarkerSize', 20);
plot(mass_outline(idx2(1:end-1),1), mass_outline(idx2(1:end-1),2), 'r.', 'MarkerSize', 20);
plot(mass_outline(idx1(1:end-1),1), mass_outline(idx1(1:end-1),2), 'rx', 'MarkerSize', 20);
plot(mass_outline(idx2(1:end-1),1), mass_outline(idx2(1:end-1),2), 'rx', 'MarkerSize', 20);
plot(mass_outline(m_idx1,1), mass_outline(m_idx1,2), 'g.', 'MarkerSize', 25);
plot(mass_outline(m_idx2,1), mass_outline(m_idx2,2), 'g.', 'MarkerSize', 25);
plot(mass_outline(m_idx1,1), mass_outline(m_idx1,2), 'gx', 'MarkerSize', 25);
plot(mass_outline(m_idx2,1), mass_outline(m_idx2,2), 'gx', 'MarkerSize', 25);

plot(mass_centroid(:,1), mass_centroid(:,2), 'bx', 'MarkerSize', 20);
% plot([mass_centroid(:,1) mass_outline(m_idx2,1)], [mass_centroid(:,2) mass_outline(m_idx2,2)], 'b:', 'LineWidth', 2);

% annotation('arrow',...
%     'Units', 'pixels',...
%     'Position', [mass_centroid - [x_min y_min] 1.5*dist*c_to_n(1,:)],...
%     'LineStyle', ':',...
%     'LineWidth', 2,...
%     'Color', 'b');
% 
% text(...
%     'Position', p1 + [0 10],...
%     'String', 'To the nipple',...
%     'VerticalAlignment', 'bottom',...
%     'HorizontalAlignment', 'right',...
%     'FontSize', 16);
% 
% text(...
%     'Position', mass_centroid + [20 10],...
%     'String', 'centroid of the mass',...
%     'VerticalAlignment', 'bottom',...
%     'HorizontalAlignment', 'right',...
%     'FontSize', 16);
% 
% text(...
%     'Position', mass_outline(1,:),...
%     'String', 'p_{1}',...
%     'VerticalAlignment', 'bottom',...
%     'HorizontalAlignment', 'left',...
%     'FontSize', 16);
% 
% text(...
%     'Position', mass_outline(idx1(2),:)+[0 5],...
%     'String', 'p_{n}',...
%     'VerticalAlignment', 'bottom',...
%     'HorizontalAlignment', 'center',...
%     'FontSize', 16);
% 
% text(...
%     'Position', mass_outline(idx2(end-1),:),...
%     'String', 'p_{2}',...
%     'VerticalAlignment', 'bottom',...
%     'HorizontalAlignment', 'left',...
%     'FontSize', 16);
% 
% text(...
%     'Position', mass_outline(m_idx2,:)-[5 0],...
%     'String', 'p_{n/2 + 1}',...
%     'VerticalAlignment', 'middle',...
%     'HorizontalAlignment', 'right',...
%     'FontSize', 16);

print('-dtiff', '-noui', '-painters', ['-f' num2str(f1)], '-r300', 'C:\isbe\thesis\figures\6\mass_border_points');
%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Outline of a breast border
borders_list = dir('C:\isbe\dev\segmentation\breast_borders\*ML*.mat');
%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Computing component weights
load C:\isbe\dev\mass_model\weights\error_map.mat

figure; hold on;
mesh(shape_points, tex_points, error_map);
xlabel('W_shape');
ylabel('W_tex');
zlabel('MGE');

[shape_i tex_i] = meshgrid(linspace(0,0.005,100), linspace(0,0.001,100));
error_map_i = interp2(shape_points, tex_points, error_map, shape_i, tex_i, 'bicubic');
figure; contour(shape_i,tex_i,error_map_i);
figure; mesh(shape_i,tex_i,error_map_i);
%print('-dtiff', '-noui', '-painters', gcf, '-r300', 'C:\isbe\thesis\figures\6\weigths_error_map.tif');

figure; mesh(shape_points, tex_points, error_map);
%print('-dtiff', '-noui', '-painters', gcf, '-r300', 'C:\isbe\thesis\figures\6\weigths_error_map_coarse.tif');

[min_er, im] = min(error_map(:));
hold on;
[sp tp] = meshgrid(shape_points, tex_points);
plot3(sp(im), tp(im), min_er, 'rx');
%%
load C:\isbe\dev\mass_model\weights\er_var
load C:\isbe\dev\mass_model\weights\er_sd
load C:\isbe\dev\mass_model\weights\er_fd
load C:\isbe\dev\mass_model\weights\er_opt3.mat

%%
% Boxplot of errors
figure('windowstyle', 'normal', 'position', [50 50 600 300], 'PaperPositionMode','auto');
hold on;
plot([mean(er_var.com.weights), mean(er_var.com.weights)], [3.75, 4.25], 'r', 'LineWidth', 1.5);
plot([mean(er_var.com.weights), mean(er_var.com.weights)], [3.75, 4.25], 'g', 'LineWidth', 1.5);
legend({'Median', 'Mean'}, 'location', 'southeast');
plot([mean(er_sd.com.weights), mean(er_sd.com.weights)], [2.75, 3.25], 'g', 'LineWidth', 1.5);
plot([mean(er_fd.com.weights), mean(er_fd.com.weights)], [1.75, 2.25], 'g', 'LineWidth', 1.5);
plot([mean(er_opt3.com.weights), mean(er_opt3.com.weights)], [0.75, 1.25], 'g', 'LineWidth', 1.5);


boxplot([er_opt3.com.weights er_fd.com.weights er_sd.com.weights er_var.com.weights],...
    'whisker', 4, 'notch', 'on', 'orientation', 'horizontal');

set(gca,'YTickLabel',{'4';'3';'2';'1'})
set(gca, 'Xlim', [0 10]);
xlabel([]);
ylabel('Method');
print('-dtiff', '-noui', '-painters', gcf, '-r300', 'C:\isbe\thesis\figures\6\weights_error_box.tif');
% 
% plot([3.5, 4], [16 16], 'g', 'LineWidth', 1.5);
% plot([3.5, 4], [15 15], 'r', 'LineWidth', 1.0);
%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Synthetic masses sampled from the model
%--------------------------------------------------------------------------
load C:\isbe\dev\mass_model\models\model_w500_50K.mat
sample_args.TargetRegion = zeros(1024);
sample_args.TargetCentre = [512 512];
sample_args.MassModel = mass_model;

for ii = 1:20
    mass_name = ['C:\isbe\thesis\figures\6\sample_mass', zerostr(ii,2), '.bmp'];
    [mass] = sample_new_mass_in_region(sample_args);
    mass_region = mass.subtract_ROI(65:end-64, 65:end-64); 
    write_im_from_colormap(mass_region, mass_name, gray(256), [0 50]);
    figure; imagesc(mass_region); axis image; colormap(gray(256));
end
%%
load C:\isbe\dev\mass_model\models\model_w500_50K.mat
norm_list = dir('C:\isbe\dev\background\images\normal1024\*.bmp');
norm = double(imread(['C:\isbe\dev\background\images\normal1024\', norm_list(2).name]));

sample_args.TargetRegion = norm;
sample_args.TargetCentre = [512 512];
sample_args.MassModel = mass_model;
[mass] = sample_new_mass_in_region(sample_args);
%write_im_from_colormap(mass_region, mass_name, gray(256), [0 50]);
figure; 
subplot(1,2,1); imagesc(mass.subtract_ROI); axis image; colormap(gray(256));
subplot(1,2,2); imagesc(mass.mass_ROI); axis image; colormap(gray(256));
%%
for ii = 1:30
    norm = double(imread(['C:\isbe\dev\background\images\normal1024\', norm_list(ii).name]));
    norm_mass = norm + mass.subtract_ROI;
    figure; 
    subplot(1,2,1); imagesc(norm); axis image; colormap(gray(256));
    subplot(1,2,2); imagesc(norm_mass); axis image; colormap(gray(256));
end
%%
mass_list = dir('C:\isbe\dev\masses1024x1024\*.mat');
for ii = 31:60
    mass_bg = u_load(['C:\isbe\dev\masses1024x1024\', mass_list(ii).name]);
    try
        mass_bg_mass = mass_bg.background_ROI + mass.subtract_ROI;
        figure; 
        subplot(1,2,1); imagesc(mass_bg.background_ROI); axis image; colormap(gray(256));
        subplot(1,2,2); imagesc(mass_bg_mass); axis image; colormap(gray(256));
    catch
        figure;
    end
end
%%
load C:\isbe\thesis\figures\6\sample_mass.mat mass
norm = double(imread(['C:\isbe\dev\background\images\normal1024\', norm_list(12).name]));
norm_mass = norm + mass.subtract_ROI;
mass_bg = u_load(['C:\isbe\dev\masses1024x1024\', mass_list(44).name]);
mass_bg_mass = mass_bg.background_ROI + mass.subtract_ROI;

mass_bg_mass = mass_bg_mass(65:end-64, 65:end-64);
norm_mass = norm_mass(65:end-64, 65:end-64);
%
c_min = min([norm_mass(:); mass_bg_mass(:)]);
c_max = max([norm_mass(:); mass_bg_mass(:)]);
%
figure; imagesc(norm_mass); axis image; colormap(gray(256)); caxis([c_min c_max]);
figure; imagesc(mass_bg_mass); axis image; colormap(gray(256));  caxis([c_min c_max]);

write_im_from_colormap(mass.subtract_ROI, 'C:\isbe\thesis\figures\6\next_ch_1.bmp', gray(256), [0 50]);
write_im_from_colormap(norm(65:end-64, 65:end-64), 'C:\isbe\thesis\figures\6\next_ch_2.bmp', gray(256), [c_min c_max]);
write_im_from_colormap(mass_bg.background_ROI(65:end-64, 65:end-64), 'C:\isbe\thesis\figures\6\next_ch_3.bmp', gray(256), [c_min c_max]);
write_im_from_colormap(norm_mass, 'C:\isbe\thesis\figures\6\next_ch_4.bmp', gray(256), [c_min c_max]);
write_im_from_colormap(mass_bg_mass, 'C:\isbe\thesis\figures\6\next_ch_5.bmp', gray(256), [c_min c_max]);
%%
load C:\isbe\dev\files\u_files.mat
[model_original_shapes mlo_id]...
    = generate_mass_AM(u_files1, 'C:\isbe\dev\mass_model\models\model_original_shapes', 'shiftOrigin', 0);
%%
[er_opt2.com er_opt2.ind] = ...
    model_errors_loo2(mass_model, u_files1, 'weights', weights_opt2);