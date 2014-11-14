%--------------------------------------------------------------------------
% ILP theory figures
%--------------------------------------------------------------------------

%Fig 1: Build 1D step-edges
step_edges1D = 2*ones(64, 64);
for offset = 1:64;
    step_edges1D(1:offset, offset) = 1;
end

coeff = 4;
level = 3;
dim = 2^level;
st = (dim+1)/2;


%Construct 1D dual-tree for each offset step-edge
[dummy dt] = dtwavexfm(step_edges1D, 4, 'near_sym_b','qshift_b');
%
%-------------------------------------------------------------------------
f1 = figure('Position', [100,100, 400, 400], 'WindowStyle', 'normal', 'Color', [1, 1, 1],'PaperPositionMode','auto');

%create axes to show the profile of the level phases
a1 = axes('Units', 'pixels', 'position', [50, 50, 300, 300]); hold on; axis([5 49 -pi pi]);

step_offset1 = 32;

%plot level 4 phase
plot(5:19, angle(dt{level+1}(coeff/2,1:15)), 'g', 'LineWidth', 2.0);
plot(20:42, angle(dt{level+1}(coeff/2,16:38)), 'g', 'LineWidth', 2.0);
plot(43:64, angle(dt{level+1}(coeff/2,39:60)), 'g', 'LineWidth', 2.0);
level4_phase = angle(dt{level+1}(coeff/2,step_offset1-4));
theta4 = text(step_offset1, level4_phase, '\theta_4');

%plot level 3 phase
plot(1:23, angle(dt{level}(coeff,1:23)), 'r', 'LineWidth', 2.0);
plot(24:35, angle(dt{level}(coeff,24:35)), 'r', 'LineWidth', 2.0);
plot(36:64, angle(dt{level}(coeff,36:64)), 'r', 'LineWidth', 2.0);
level3_phase = angle(dt{level}(coeff,step_offset1));
theta3 = text(step_offset1, level3_phase, '\theta_3');

title('DT-CWT coefficient phase (\theta) vs step offset (x)');
xlabel('Offset of step from coefficient  (x)')
ylabel('DT-CWT phase (\theta)');

set(a1, 'Xtick', st+(0:7)*dim, 'Xticklabel', (-3:4)*dim, 'Ytick', [-pi 0 pi], 'Yticklabel', {'-Pi', '0', 'Pi'});

%print_eps('K:\isbe\thesis\figures\7\ilp_theory_1.eps');
print('-dtiff', '-noui', '-painters', f1, '-r300', '\isbe\thesis\figures\7\ilp_theory_1.tif');
%
%-------------------------------------------------------------------------
f1 = figure('Position', [500,500, 400, 400], 'WindowStyle', 'normal', 'Color', [1, 1, 1],'PaperPositionMode','auto');

%create axes to show the profile of the level phases
a1 = axes('Units', 'pixels', 'position', [50, 50, 300, 300]); hold on; axis([5 49 -pi pi]);

step_offset1 = 32;

%plot level 4 phase
plot(5:11, angle(dt{level+1}(coeff/2,1:7).^2), 'g', 'LineWidth', 2.0);
plot(12:25, angle(dt{level+1}(coeff/2,8:21).^2), 'g', 'LineWidth', 2.0);
plot(26:36, angle(dt{level+1}(coeff/2,22:32).^2), 'g', 'LineWidth', 2.0);
plot(37:64, angle(dt{level+1}(coeff/2,33:60).^2), 'g', 'LineWidth', 2.0);
level4_phase = angle(dt{level+1}(coeff/2,step_offset1-4).^2);
theta4 = text(step_offset1, level4_phase, '\theta_4');

%plot level 3 phase
plot(1:23, angle(dt{level}(coeff,1:23)), 'r', 'LineWidth', 2.0);
plot(24:35, angle(dt{level}(coeff,24:35)), 'r', 'LineWidth', 2.0);
plot(36:64, angle(dt{level}(coeff,36:64)), 'r', 'LineWidth', 2.0);
level3_phase = angle(dt{level}(coeff,step_offset1));
theta3 = text(step_offset1, level3_phase, '\theta_3');

title('DT-CWT coefficient phase (\theta) vs step offset (x)');
xlabel('Offset of step from coefficient  (x)')
ylabel('DT-CWT phase (\theta)');

set(a1, 'Xtick', st+(0:7)*dim, 'Xticklabel', (-3:4)*dim, 'Ytick', [-pi 0 pi], 'Yticklabel', {'-Pi', '0', 'Pi'});
%print_eps('K:\isbe\thesis\figures\7\ilp_theory_2.eps');
print('-dtiff', '-noui', '-painters', f1, '-r300', '\isbe\thesis\figures\7\ilp_theory_2.tif');

%%
%--------------------------------------------------------------------------
% Original ICP ranges over the circle
[...
    -pi/4.49        pi/4.49         (180/pi)*(-pi/4.49)     (180/pi)*(pi/4.49);...
    pi/4-pi/8.89    pi/4+pi/8.89    (180/pi)*(pi/4-pi/8.89) (180/pi)*(pi/4+pi/8.89);...
    pi/2-pi/4.49    pi/2+pi/4.49    (180/pi)*(pi/2-pi/4.49) (180/pi)*(pi/2+pi/4.49);...
    pi/2-pi/4.49    pi/2+pi/4.49    (180/pi)*(pi/2-pi/4.49) (180/pi)*(pi/2+pi/4.49);...
    3*pi/4-pi/8.89  3*pi/4+pi/8.89  (180/pi)*(3*pi/4-pi/8.89) (180/pi)*(3*pi/4+pi/8.89);...
    -pi/4.49        pi/4.49         (180/pi)*(-pi/4.49)     (180/pi)*(pi/4.49);...
]
%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Detecting linear structures figures
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Create ILP/ICP images for real mammogram plots, can we show a normal
% patch and a radial patch?
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%Find streaky patch from ILP report

%try backgrounds, 006, 008 and 004 for radial texture
bg = double(imread('C:\isbe\dev\background\images\mass_2\mass008.bmp'));
%bg = bg(142:1100, :);
dt_bg = dtwavexfm2(bg, 7);

[ilp_bg icp_bg] = mb_dual_tree_transform(dt_bg);
dt_bg2 = dt_bg;
%
for level = 2:5
    
    [max_ilp{level}] = max(ilp_bg{level}, [], 3);
    %[max_icp{level}] = max(icp_bg{level}, [], 3);
    [max_icp{level}] = compute_max_icp(icp_bg{level});
    
    ilp_mag = histeq(abs(max_ilp{level}) / max(abs(max_ilp{level}(:))), 1 ./ (1:256));
    icp_mag = histeq(abs(max_icp{level}) / max(abs(max_icp{level}(:))), 1 ./ (1:256));
    temp = icp_mag.*exp(i*angle(max_icp{level}));
    figure; 
    image(complex2rgb(ilp_mag.*exp(i*2*angle(max_icp{level})))); axis image; hold on;
    title(['ICP level ', num2str(level)]);
    %quiver(real(max_icp{level}), -imag(max_icp{level}), 'w');
    %figure;
    %quiver(real(temp), -imag(temp)); axis image;
    figure; 
    image(complex2rgb(ilp_mag.*exp(i*angle(max_ilp{level})))); axis image;
    hold on; title(['ILP level ', num2str(level)]);
    
    [bw] = icp_hysterisis(max_icp{level}, [cutoffs(level,80) cutoffs(level,95)]);
    
%     bw = ilp_hysterisis(max_ilp{level}, [0.075 0.01], [-inf 0]);
    icp_mag(~bw) = 0;
    figure; title(['ILP level ', num2str(level)]);
    image(complex2rgb(double(icp_mag).*exp(2*i*angle(max_icp{level})))); axis image;
%     figure; 
%     image(complex2rgb(double(ilp_mag).*exp(i*angle(max_ilp{level})))); axis image;
%     title(['ILP level ', num2str(level)]);
    
    for ori = 1:6
        temp = dt_bg2{level}(:,:,ori);
        temp(~bw) = 0;
        dt_bg2{level}(:,:,ori) = temp;
    end
    %image(complex2rgb(double(ilp_mag).*exp(i*2*angle(max_icp{level}))));
    %axis image;
   % quiver(real(max_ilp{level}), -imag(max_ilp{level}), 'w');
   
%     imwrite(complex2rgb(ilp_mag.*exp(i*angle(max_ilp{level}))),...
%         ['C:\isbe\thesis\figures\7\ilp_bg_mp', num2str(level), '.bmp']);
%     imwrite(complex2rgb(ilp_mag.*exp(i*2*angle(max_icp{level}))),...
%         ['C:\isbe\thesis\figures\7\icp_bg_mp', num2str(level), '.bmp']);
end
%
dt_bg2{6} = zeros(size(dt_bg2{6}));
dt_bg2{7} = zeros(size(dt_bg2{7}));
line_image = dtwaveifm2(dt_bg2, 2);
figure; imagesc(line_image); axis image; colormap(gray(256));
%%
write_im_from_colormap(bg, 'C:\isbe\thesis\figures\7\mass008.bmp', gray(256));
%%
%--------------------------------------------------------------------------
% Streaky background
%--------------------------------------------------------------------------
streaky = rot90(double(imread('C:\isbe\dev\background\images\normal512\o04_011RML_1024_3309_2896.bmp')));
streaky_dt = dtwavexfm2(streaky, 7);

[streaky_ilp streaky_icp] = mb_dual_tree_transform(streaky_dt);
dt_bg2 = streaky_dt;
%
for level = 2:5
    
    [max_ilp{level}] = max(streaky_ilp{level}, [], 3);
    [max_icp{level}] = max(streaky_icp{level}, [], 3);
    
    ilp_mag = histeq(abs(max_ilp{level}) / max(abs(max_ilp{level}(:))), 1 ./ (1:256));
    icp_mag = histeq(abs(max_icp{level}) / max(abs(max_icp{level}(:))), 1 ./ (1:256));
    temp = icp_mag.*exp(i*angle(max_icp{level}));
    figure; title(['ICP level ', num2str(level)]);
    image(complex2rgb(ilp_mag.*exp(i*2*angle(max_icp{level})))); axis image; hold on;
    %quiver(real(max_icp{level}), -imag(max_icp{level}), 'w');
    %figure;
    %quiver(real(temp), -imag(temp)); axis image;
    figure; title(['ILP level ', num2str(level)]);
    image(complex2rgb(ilp_mag.*exp(i*angle(max_ilp{level})))); axis image;
    hold on;
   % quiver(real(max_ilp{level}), -imag(max_ilp{level}), 'w');
   
   [bw] = icp_hysterisis(max_icp{level}, [cutoffs(level,80) cutoffs(level,95)]);
   imwrite(double(bw),...
        ['C:\isbe\thesis\figures\7\bw_streaky_mp', num2str(level), '.bmp']);
%     bw = ilp_hysterisis(max_ilp{level}, [0.07 0.09], [-inf -inf]);
    icp_mag(~bw) = 0;
    figure; title(['ILP level ', num2str(level)]);
    image(complex2rgb(double(icp_mag).*exp(2*i*angle(max_icp{level})))); axis image;
    
    for ori = 1:6
        temp = dt_bg2{level}(:,:,ori);
        temp(~bw) = 0;
        dt_bg2{level}(:,:,ori) = temp;
    end
    
    %image(complex2rgb(double(ilp_mag).*exp(i*2*angle(max_icp{level})))); axis image;
    
%     imwrite(complex2rgb(ilp_mag.*exp(i*angle(max_ilp{level}))),...
%         ['C:\isbe\thesis\figures\7\ilp_streaky_mp', num2str(level), '.bmp']);
%     imwrite(complex2rgb(ilp_mag.*exp(i*2*angle(max_icp{level}))),...
%         ['C:\isbe\thesis\figures\7\icp_streaky_mp', num2str(level), '.bmp']);
end


dt_bg2{6} = zeros(size(dt_bg2{6}));
dt_bg2{7} = zeros(size(dt_bg2{7}));
line_image = dtwaveifm2(dt_bg2, 2);
figure; imagesc(line_image); axis image; colormap(gray(256));

%%
streaky_ilp_abs = zeros(size(streaky));
for level = 2:5
    si = imresize(abs(max_ilp{level}), size(streaky_ilp_abs), 'bilinear');
    si = si / max(si(:));
    
    figure;
    subplot(1,2,1); imagesc(streaky_ilp_abs); axis image;
    subplot(1,2,2); imagesc(si); axis image;
    streaky_ilp_abs = max(streaky_ilp_abs, si);
end
%%
%--------------------------------------------------------------------------
% Find magnitude thresholds for ICP
%--------------------------------------------------------------------------
%%
hist_counts = zeros(5,200);
data_type = 'normal512';

im_list = dir([mberksroot, 'background/images/' data_type , '/*.bmp']);
bin_centres = linspace(0,10,200);

for ii = 1:length(im_list);

    bg = double(imread([mberksroot,...
        'background/images/', data_type, '/', im_list(ii).name]));
    
    dual_tree = dtwavexfm2(bg, 7);
    
    [ilp icp] = mb_dual_tree_transform(dual_tree);
    %
    for level = 1:5

        %[max_ilp] = max(ilp{level}, [], 3);
        [max_icp] = compute_max_icp(icp{level});
        counts = hist(abs(max_icp(:)), bin_centres);
        hist_counts(level,:) = hist_counts(level,:) + counts;
    end
    clear ilp icp dual_tree bg
end
%
cutoffs = zeros(5,100);
for level = 1:5
    cum_sum = 100*cumsum(hist_counts(level,:)) / sum(hist_counts(level,:));
    for ii = 1:100
        cutoffs(level,ii) = bin_centres(find(cum_sum >= ii, 1));
    end
end
save C:\isbe\dev\background\cls\thresholds cutoffs hist_counts bin_centres
%%
[template_orientation_map] =  mb_dt_orientation_map('Image', streaky);
n_bins = 120;
figure; weighted_complex_rose(template_orientation_map(:), n_bins);
%%
[template_orientation_map] =  mb_dt_orientation_map('Image', bg);
n_bins = 120;
figure; weighted_complex_rose(template_orientation_map(:), n_bins);
%%
[template_ilp_map] =  mb_dt_ilp_map('Image', streaky);
n_bins = 120;
figure; weighted_complex_rose(template_ilp_map(:), n_bins);
%
[template_ilp_map] =  mb_dt_ilp_map('Image', bg);
n_bins = 120;
figure; weighted_complex_rose(template_ilp_map(:), n_bins);
%%
% What does an image of only the maximal coefficients look like?
streaky = rot90(double(imread('C:\isbe\dev\background\images\normal512\o04_011RML_1024_3309_2896.bmp')));
streaky_dt = dtwavexfm2(streaky, 7);
for level = 1:7
    [dummy max_idx] = max(streaky_dt{level}, [], 3);
    
    for ori = 1:6
        band = streaky_dt{level}(:,:,ori);
        band(max_idx ~= ori) = 0;
        streaky_dt{level}(:,:,ori) = band;
    end
end
streaky_max = dtwaveifm2(streaky_dt);
figure;
subplot(1,2,1); imagesc(streaky); axis image; colormap(gray(256));
subplot(1,2,2); imagesc(streaky_max); axis image; colormap(gray(256));
write_im_from_colormap(streaky_max, 'C:\isbe\thesis\figures\7\streaky_max.bmp', gray(256));
%%
bg = double(imread('C:\isbe\dev\background\images\mass_2\mass008.bmp'));
bg = bg(71:550, :);
bg_dt = dtwavexfm2(bg, 7);
for level = 1:7
    [dummy max_idx] = max(bg_dt{level}, [], 3);
    
    for ori = 1:6
        band = bg_dt{level}(:,:,ori);
        band(max_idx ~= ori) = 0;
        bg_dt{level}(:,:,ori) = band;
    end
end
bg_max = dtwaveifm2(bg_dt);
figure;
subplot(1,2,1); imagesc(bg); axis image; colormap(gray(256));
subplot(1,2,2); imagesc(bg_max); axis image; colormap(gray(256));
write_im_from_colormap(bg_max, 'C:\isbe\thesis\figures\7\bg_max.bmp', gray(256));
%%

%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% DT Synthesis figures

%Models built by submitting band_cluster_dual_tree_all to hydra for levels 1
%to 4

% Testing the models:
data_type = 'normal512';
dt_list = dir([mberksroot, 'background/dual_tree/' data_type , '/*dual_tree*']);


display(['--Texture synthesis script started: ' datestr(now)]);

for ii = 1:1
    
    if ii > length(dt_list); break; end;
    
    dt_args.TargetDualTree = u_load([mberksroot,...
        'background/dual_tree/', data_type, '/', dt_list(ii).name]);
    
    [rows cols] = size(dt_args.TargetDualTree{1}(:,:,1));
    row_centre = round(rows / 2);
    col_centre = round(cols / 2);

    % Make a the biggest circular mask
    m = min([rows cols]);
    rad = floor((m - 128) / 2);
    [x y] = meshgrid(1:cols, 1:rows);
    dt_args.FilledImage = (x - col_centre).^2 + (y - row_centre).^2 > rad.^2;
    clear x y m rad row_centre col_centre;

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dt_args.ModelDir = [mberksroot, 'background/models/dual_tree/'];
    dt_args.ModelName = 'mass_2_k10_w1_3_w2_3';
    dt_args.SynthesisMode = 'patch-wise';
    dt_args.MovieFilename = [mberksroot, 'background/syn/dual_tree/',...
        data_type, '_', dt_args.ModelName, '_' zerostr(ii, 3), '.gif'];
    dt_args.SaveFile = [mberksroot, 'background/syn/dual_tree/',...
        data_type, '_', dt_args.ModelName, '_' zerostr(ii, 3)];

    %[synthesised_image, dual_tree] = mb_gmm_dual_tree_synthesis_final(dt_args);
    [synthesised_image, dual_tree] = mb_gmm_dual_tree_synthesis_oris2(dt_args);
    
    %clear dt_args;
end

display(['--Texture synthesis script finished: ' datestr(now)]);

%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Correlations within the dual-tree
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%%
%Get indices of points to sample
clear;
DualTreeDir = 'C:\isbe\dev\background\dual_tree\normal512\';
WindowSize = 11;
MaxMemory = 256; 
overlap = 10;
figure;
for Level = 1:4;

    % Set constant
    size_of_data_element = 8; % assume each data element is a double, hence 8 bytes per element.

    % first get a directory listing and other info about the image directory
    DualTreeFiles = dir([DualTreeDir,'*dual_tree*']);
    number_of_trees = length(DualTreeFiles);

    % work out how much memory (in bytes) each window will occupy
    mem_per_sample = WindowSize^2 * size_of_data_element;

    % work out how many points we can fit in MaxMemory
    num_points = floor((MaxMemory * 1024 * 1024) /  mem_per_sample);

    % work out how many points per image this is
    num_points_per_tree = floor(num_points / number_of_trees);

    %display(['Next Data Function: num_points_per_image = ' num2str(num_points_per_tree)])

    all_indices = repmat(NaN, number_of_trees, num_points_per_tree);

    for ii = 1:number_of_trees

        % open the i-th image file
        dual_tree = u_load([DualTreeDir, DualTreeFiles(ii).name]);

        % create a matrix of zeros of the same size as the image
        sampled = zeros(size(dual_tree{Level}(:,:,1)));

        % set the edges to ones so that we don't sample windows that lie off
        % the edge of the window
        sampled(1:overlap, :) = 1; % top portion
        sampled(end-overlap+1:end, :) = 1; % bottom portion
        sampled(:, 1:overlap) = 1; % left side
        sampled(:,end-overlap+1 : end) = 1;%#ok % right side    

        % sample from the image
        unsampled_indices = find(~sampled);
        points_in_this_tree = min(num_points_per_tree, length(unsampled_indices));

        rp = randperm(length(unsampled_indices));    
        all_indices(ii, 1:points_in_this_tree) = ...
            unsampled_indices(rp(1:points_in_this_tree));

    end


    data_args.LoadIdx = 0;
    data_args.Idx = all_indices; clear all_indices;
    data_args.DualTreeDir = DualTreeDir;
    data_args.Level = Level;
    data_args.WindowSize = WindowSize;
    data_args.WindowSize2 = 0;
    data_args.Standardise = 0;
    data_args.Real = 1;
    data_args.Orientation = 1; %
    [data] = mb_get_dual_tree_subband_data(data_args);
    
%     cov_mat{Level} = cov(data1); %#ok
%     sd_mat = sqrt(diag(cov_mat{Level})) * ones(1, size(data1,2));
%     corr_mat{Level} = cov_mat{Level} ./ (sd_mat .* sd_mat'); %#ok
    [corr_mat{Level},p_mat{Level}] = corr(data);
    corr_mat{Level}(logical(eye(121))) = 0;
    clear data;
% %%    
% figure;
% for Level = 1:4
    %subplot(2,2,Level); imagesc(reshape(abs(corr_mat{Level}(61,:)),11,11)); axis image;
    subplot(2,2,Level); imagesc(reshape(p_mat{Level}(61,:),11,11)); axis image;
    write_im_from_colormap(kron(reshape(abs(corr_mat{Level}(61,:)),11,11), ones(16)),...
        ['C:\isbe\thesis\figures\7\corr_map_local_dt_', num2str(Level), '.bmp'], ...
        jet(256));
    write_im_from_colormap(kron(reshape(p_mat{Level}(61,:)>.05,11,11), ones(16)),...
        ['C:\isbe\thesis\figures\7\p_map_local_dt_', num2str(Level), '.bmp'], ...
        jet(256));
    display(['Maximum local correlation for level ', num2str(Level), ' = ', num2str(max(abs(corr_mat{Level}(61,:))))]);
end

%%
clear;
DualTreeDir = 'C:\isbe\dev\background\dual_tree\normal512\';
WindowSize = 11;
MaxMemory = 32; 
overlap = 10;
figure;
for Level = 1:4;

    % Set constant
    size_of_data_element = 8; % assume each data element is a double, hence 8 bytes per element.

    % first get a directory listing and other info about the image directory
    DualTreeFiles = dir([DualTreeDir,'*dual_tree*']);
    number_of_trees = length(DualTreeFiles);

    % work out how much memory (in bytes) each window will occupy
    mem_per_sample = 12 * size_of_data_element;

    % work out how many points we can fit in MaxMemory
    num_points = floor((MaxMemory * 1024 * 1024) /  mem_per_sample);

    % work out how many points per image this is
    num_points_per_tree = floor(num_points / number_of_trees);

    %display(['Next Data Function: num_points_per_image = ' num2str(num_points_per_tree)])

    all_indices = repmat(NaN, number_of_trees, num_points_per_tree);

    for ii = 1:number_of_trees

        % open the i-th image file
        dual_tree = u_load([DualTreeDir, DualTreeFiles(ii).name]);

        % create a matrix of zeros of the same size as the image
        sampled = zeros(size(dual_tree{Level}(:,:,1)));

        % set the edges to ones so that we don't sample windows that lie off
        % the edge of the window
        sampled(1:overlap, :) = 1; % top portion
        sampled(end-overlap+1:end, :) = 1; % bottom portion
        sampled(:, 1:overlap) = 1; % left side
        sampled(:,end-overlap+1 : end) = 1;%#ok % right side    

        % sample from the image
        unsampled_indices = find(~sampled);
        points_in_this_tree = min(num_points_per_tree, length(unsampled_indices));

        rp = randperm(length(unsampled_indices));    
        all_indices(ii, 1:points_in_this_tree) = ...
            unsampled_indices(rp(1:points_in_this_tree));

    end


    data_args.LoadIdx = 0;
    data_args.Idx = all_indices; clear all_indices;
    data_args.DualTreeDir = DualTreeDir;
    data_args.Level = Level;
    data_args.WindowSize = 1;
    data_args.WindowSize2 = 0;
    data_args.Standardise = 0;
    
    for ori = 1:6
        data_args.Real = 1;
        data_args.Orientation = ori; %#ok
        data(:,ori) = mb_get_dual_tree_subband_data(data_args);
        
        data_args.Real = 0;
        data_args.Orientation = ori; %#ok
        data(:,ori+6) = mb_get_dual_tree_subband_data(data_args);
    end
    
%     cov_mat{Level} = cov(data); %#ok
%     sd_mat = sqrt(diag(cov_mat{Level})) * ones(1, size(data,2));
%     corr_mat{Level} = cov_mat{Level} ./ (sd_mat .* sd_mat'); %#ok
    [corr_mat{Level},p_mat{Level}] = corr(data);
    corr_mat{Level}(logical(eye(12))) = 0;
    clear data;
    
% %%    
% figure;
% for Level = 1:4
    %subplot(2,2,Level); imagesc(abs(corr_mat{Level})); axis image;
    subplot(2,2,Level); imagesc(p_mat{Level} <.05); axis image;
    
    write_im_from_colormap(kron(abs(corr_mat{Level}), ones(16)),...
        ['C:\isbe\thesis\figures\7\corr_map_oris_dt_', num2str(Level), '.bmp'], ...
        jet(256));
    
    display(['Maximum orientation correlation for level ', num2str(Level), ' = ', num2str(max(abs(corr_mat{Level}(:))))]);
end
%%
clear;
DualTreeDir = 'C:\isbe\dev\background\dual_tree\normal512\';
MaxMemory = 32; 
overlap = 10;    

% Set constant
size_of_data_element = 8; % assume each data element is a double, hence 8 bytes per element.

% first get a directory listing and other info about the image directory
DualTreeFiles = dir([DualTreeDir,'*dual_tree*']);
number_of_trees = length(DualTreeFiles);

% work out how much memory (in bytes) each window will occupy
mem_per_sample = 60 * size_of_data_element;

% work out how many points we can fit in MaxMemory
num_points = floor((MaxMemory * 1024 * 1024) /  mem_per_sample);

% work out how many points per image this is
num_points_per_tree = floor(num_points / number_of_trees);

%display(['Next Data Function: num_points_per_image = ' num2str(num_points_per_tree)])

data = repmat(NaN, 12*number_of_trees*num_points_per_tree, 5);
curr_idx = 1;
for ii = 1:number_of_trees

    % open the i-th image file
    dual_tree = u_load([DualTreeDir, DualTreeFiles(ii).name]);

    % create a matrix of zeros of the same size as the image
    sampled = zeros(size(dual_tree{1}(:,:,1)));

    % set the edges to ones so that we don't sample windows that lie off
    % the edge of the window
    sampled(1:overlap, :) = 1; % top portion
    sampled(end-overlap+1:end, :) = 1; % bottom portion
    sampled(:, 1:overlap) = 1; % left side
    sampled(:,end-overlap+1 : end) = 1;%#ok % right side    

    % sample from the image
    [r1_pts c1_pts] = find(~sampled);
    points_in_this_tree = min(num_points_per_tree, length(r1_pts));

    rp = randperm(length(r1_pts));    
    r1_pts = r1_pts(rp(1:points_in_this_tree));
    c1_pts = c1_pts(rp(1:points_in_this_tree));

    for ori = 1:6
        for level = 1:5 


            r_pts = ceil(r1_pts/2^(level-1));
            c_pts = ceil(c1_pts/2^(level-1));


            idx = sub2ind(size(dual_tree{level}(:,:,1)), r_pts, c_pts);


            subband = dual_tree{level}(:,:,ori);
            data(curr_idx:curr_idx+points_in_this_tree-1,level) = real(subband(idx));             
            data(curr_idx+points_in_this_tree:curr_idx+2*points_in_this_tree-1,level) = imag(subband(idx));
        end
        curr_idx = curr_idx + 2*points_in_this_tree;
    end
end

data(curr_idx:end,:) = [];
data(any(isnan(data),2),:) = [];

cov_mat = cov(data); %#ok

sd_mat = sqrt(diag(cov_mat)) * ones(1, size(cov_mat,1));
corr_mat = cov_mat ./ (sd_mat .* sd_mat'); %#ok
[corr_mat,p_mat] = corr(data);
corr_mat(logical(eye(5))) = 0;%#ok
clear data*;

figure; imagesc(abs(corr_mat)); axis image;
figure; imagesc(p_mat < .05); axis image;

write_im_from_colormap(kron(abs(corr_mat), ones(32)),...
    ['C:\isbe\thesis\figures\7\corr_map_levels_dt.bmp'], ...
    jet(256));

display(['Maximum level correlation = ', num2str(max(abs(corr_mat(:))))]);

%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Synthesising from the models
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%%
data_type = 'normal512';
dt_list = dir([mberksroot, 'background/dual_tree/' data_type , '/*dual_tree*']);


display(['--Texture synthesis script started: ' datestr(now)]);

for ii = [12 13 21 43]%1:length(dt_list);
    
    dt_args.TargetDualTree = u_load([mberksroot,...
        'background/dual_tree/', data_type, '/', dt_list(ii).name]);
    
    [rows cols] = size(dt_args.TargetDualTree{1}(:,:,1));
    row_centre = round(rows / 2);
    col_centre = round(cols / 2);

    % Make a the biggest circular mask
    m = min([rows cols]);
    rad = floor((m - 64) / 2);
    [x y] = meshgrid(1:cols, 1:rows);
    dt_args.FilledImage = (x - col_centre).^2 + (y - row_centre).^2 > rad.^2;%
    clear x y m rad row_centre col_centre;

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dt_args.ModelDir = [mberksroot, 'background/models/dual_tree/'];
    dt_args.ModelName = 'normal512_k20_rep_10_w1_3_w2_1';
    dt_args.SynthesisMode = 'pixel-wise';
    dt_args.MovieFilename = [];%'test_movie.gif';
    dt_args.SaveFile = [];%[mberksroot, 'background/syn/dual_tree/',...
        %dt_args.SynthesisMode, '/',data_type, '_', dt_args.ModelName, '_a' zerostr(ii, 3)];

    [synthesised_image, dual_tree] = mb_gmm_dual_tree_synthesis_final(dt_args);
    syn_name = ['C:\isbe\thesis\figures\8\bgs\syn_norm', zerostr(ii,3), '.bmp'];
    %write_im_from_colormap(synthesised_image, syn_name, colormap(gray(256)));
    figure; imagesc(synthesised_image); axis image; colormap(gray(256));
    %clear dt_args;
end

display(['--Texture synthesis script finished: ' datestr(now)]);
%%
data_type = 'mass_2';
dt_list = dir([mberksroot, 'background/dual_tree/' data_type , '/*dual_tree*']);


display(['--Texture synthesis script started: ' datestr(now)]);

for ii = 27%1:length(dt_list);
    
    dt_args.TargetDualTree = u_load([mberksroot,...
        'background/dual_tree/', data_type, '/', dt_list(ii).name]);
    
    [rows cols] = size(dt_args.TargetDualTree{1}(:,:,1));
    row_centre = round(rows / 2);
    col_centre = round(cols / 2);

    % Make a the biggest circular mask
    m = min([rows cols]);
    rad = floor((m - 128) / 2);
    [x y] = meshgrid(1:cols, 1:rows);
    dt_args.FilledImage = false(rows,cols);%(x - col_centre).^2 + (y - row_centre).^2 > 40.^2;%
    clear x y m rad row_centre col_centre;

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dt_args.ModelDir = [mberksroot, 'background/models/dual_tree/'];
    dt_args.ModelName = 'mass_2_k10_w1_3_w2_1';
    dt_args.SynthesisMode = 'pixel-wise';
    dt_args.MovieFilename = [];%'test_movie.gif';
    dt_args.SaveFile = [mberksroot, 'background/syn/dual_tree/',...
        dt_args.SynthesisMode, '/',data_type, '_', dt_args.ModelName, '_a' zerostr(ii, 3)];

    [synthesised_image, dual_tree] = mb_gmm_dual_tree_synthesis_frog(dt_args);

    %clear dt_args;
end

display(['--Texture synthesis script finished: ' datestr(now)]);
%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Exploring the correlations in synthetic data
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Generating CLS maps of each region in the data
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
mkdir C:\isbe\dev\background\cls\mass
load C:\isbe\dev\background\cls\thresholds
mam_list = dir('C:\isbe\dev\background\images\normal512\*.bmp');
%mam_list = dir('C:\isbe\dev\masses1024x1024\*.mat');
%
for ii = 1:length(mam_list)
    mam_region = double(imread(['C:\isbe\dev\background\images\normal512\', mam_list(ii).name]));
    %mass = u_load(['C:\isbe\dev\masses1024x1024\', mam_list(ii).name]);
    %mam_region = double(mass.mass_ROI);
    [cls] = mb_dual_tree_cls_detect(mam_region, 5, cutoffs, [90 10], [-inf pi/5], 0);
    cls_name = ['C:\isbe\dev\background\cls\normal512\', mam_list(ii).name(1:end-4), '_cls.mat'];
    save(cls_name, 'cls');
    display(['processing: ', num2str(ii), ' of ', num2str(length(mam_list))]);
    
%     figure;
%     imagesc(mam_region); axis image; colormap(gray(256)); hold on;
%     
%     colors = 'rrrgb';
%     for jj = 3:5
%         [r c] = find(bwmorph(imresize(cls.map{jj}, size(mam_region)), 'thin', inf));
%         plot(c,r,[colors(jj), '.'], 'MarkerSize', 4);
%     end 
%     figure;
%     subplot(2,2,1); imagesc(mam_region); axis image; colormap(gray(256));
%     subplot(2,2,2); imagesc(cls.map{3}); axis image; colormap(gray(256));
%     subplot(2,2,3); imagesc(cls.map{4}); axis image; colormap(gray(256));
%     subplot(2,2,4); imagesc(cls.map{5}); axis image; colormap(gray(256));
    clear mam_region cls
end
%%
for ii = 28:29;% 24 34 36 40 51 52 53]
    
    %for each image
    mass = u_load(['C:\isbe\dev\masses1024x1024\', mam_list(ii).name]);
    mam_region = double(mass.mass_ROI);
    
    [gp p_sizes] = buildGpyr(mam_region, 6);
    gp = mb_change_pyramid_form(gp, p_sizes, 'g');
    
    CLS_results = cell(3,1);
    gf_args.NumAngles = 36;
    for jj = 3:5
        %Do CLS detection
        [CLS_results{jj-2}] = mb_cls_selection('ImageIn', gp{jj}, 'GaborFilterArgs', gf_args);
    end
    
%     save(['C:\isbe\dev\background\misc\mass_cls_', zerostr(ii,3)], 'CLS_results');
    
    figure('Name', ['mass', zerostr(ii,3)]); imagesc(gp{2});
    colormap(gray(256)); axis image; hold on;
    colors = 'rgb';
    for jj = 1:3
        [r c] = find(bwmorph(imresize(CLS_results{jj}.CLS, size(gp{2})), 'thin', inf));
        plot(c,r,[colors(jj), '.'], 'MarkerSize', 4);
    end
    
%     cls_mr = CLS_results{1}.CLS | ...
%         imresize(CLS_results{2}.CLS, size(CLS_results{1}.CLS)) | ...
%         imresize(CLS_results{3}.CLS, size(CLS_results{1}.CLS));
%     [r c] = find(cls_mr);    
%     figure('Name', ['mass', zerostr(ii,3)]); imagesc(gp{2});
%     colormap(gray(256)); axis image; hold on;
%     plot(c,r,'r.', 'MarkerSize', 4);
    
    clear i1 r c CLS_results cls_mr;
end
%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Synthesise using CLS models
%-------------------------------------------------------------------------
data_type = 'normal512';
load C:\isbe\dev\files\u_files.mat
cls_list = dir('C:\isbe\dev\background\cls\mass\*.mat');
dt_list = dir([mberksroot, 'background/dual_tree/' data_type , '/*dual_tree*']);

cls = u_load(['C:\isbe\dev\background\cls\mass\', cls_list(idx_u1(77)).name]);
cls.map(1) = [];
display(['--Texture synthesis script started: ' datestr(now)]);

for ii = 44%[12 13 21 43]%27%1:length(dt_list);
    
    dt_args.TargetDualTree = u_load([mberksroot,...
        'background/dual_tree/', data_type, '/', dt_list(ii).name]);
    
    [rows cols] = size(dt_args.TargetDualTree{1}(:,:,1));
    row_centre = round(rows / 2);
    col_centre = round(cols / 2);

    % Make a the biggest circular mask
    m = min([rows cols]);
    rad = floor((m - 64) / 2);
    [x y] = meshgrid(1:cols, 1:rows);
    dt_args.FilledImage = (x - col_centre).^2 + (y - row_centre).^2 > rad.^2;%
    %dt_args.FilledImage([1:32, end-31:end], :) = 1;
    %dt_args.FilledImage(:, [1:32, end-31:end]) = 1;
    clear x y m rad row_centre col_centre;

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dt_args.ModelDir = [mberksroot, 'background/models/dual_tree/'];
    dt_args.ModelName = 'normal512_k20_w1_3_w2_1';
    dt_args.SynthesisMode = 'pixel';
    dt_args.MovieFilename = [];%'test_movie.gif';
    dt_args.ClsMap = cls.map;
    dt_args.SaveFile = [mberksroot, 'background/syn/dual_tree/',...
        dt_args.SynthesisMode, '/',data_type, '_', dt_args.ModelName, '_a' zerostr(ii, 3)];

    [synthesised_image, dual_tree] = mb_gmm_dual_tree_synthesis_cls(dt_args);
    syn_name = ['C:\isbe\thesis\figures\8\syn_norm_cls', zerostr(ii,3), '.bmp'];
    write_im_from_colormap(synthesised_image, syn_name, colormap(gray(256)));
    %clear dt_args;
end

display(['--Texture synthesis script finished: ' datestr(now)]);
%%
%%
norm_list = dir('C:\isbe\dev\background\images\normal1024\*.bmp');
load C:\isbe\dev\mass_model\models\model_w500_50K.mat
mass_dir = 'C:\isbe\dev\masses1024x1024';
mkdir C:\isbe\thesis\figures\8\bgs
for nn = 16%1:length(norm_list)
    for mm = 47%[14 26 27 77]
        try
            target_region = double(imread(['C:\isbe\dev\background\images\normal1024\', norm_list(nn).name]));
            [mass] = ...
                mb_synthesise_mass(target_region, [512 512], mass_model,...
                mass_dir, model_id.mass_files, mm, []);
            mass_name = ['C:\isbe\thesis\figures\8\bgs\norm', zerostr(nn,3), 'mass', zerostr(mm,3), '.bmp'];
            write_im_from_colormap(mass, mass_name, colormap(gray(256)));
        catch
            display('Skipped');
        end
    end
%     figure;
%     subplot(1,2,1); imagesc(target_region); axis image; colormap(gray(256));
%     subplot(1,2,2); imagesc(mass.mass_ROI); axis image; colormap(gray(256));
    
end
%%
load C:\isbe\dev\files\bg_files.mat

for ii = 27%1:length(bg_files)
    mass = u_load(['C:\isbe\dev\masses1024x1024\', bg_files(ii).name]);
    figure; imagesc(mass.background_ROI); axis image; colormap(gray(256));
    clear mass;
end
%%
mass_list = u_load('C:\isbe\dev\files\u_files.mat');
norm_list = dir('C:\isbe\dev\background\images\normal1024\*.bmp');
%%
target_region = rot90(double(imread(['C:\isbe\dev\background\images\normal1024\', norm_list(16).name])),-1);
mass = u_load(['C:\isbe\dev\masses1024x1024\', mass_list(47).name]);
[mass_bg64] = ...
        mb_synthesise_mass(target_region, [512 512], mass_model,...
        mass_dir, model_id.mass_files, mm, []);
[mass_bg90] = ...
        mb_synthesise_mass(target_region, [512 512], mass_model,...
        mass_dir, model_id.mass_files, mm, 0);

%
norm_bg = target_region(129:end-128, 129:end-128);
mass_bg = mass.background_ROI(129:end-128, 129:end-128);
mass_bg64 = mass_bg64(129:end-128, 129:end-128);
mass_bg90 = mass_bg90(129:end-128, 129:end-128);
%%
figure; 
subplot(2,2,1); imagesc(norm_bg); axis image; colormap(gray(256));
subplot(2,2,2); imagesc(mass_bg); axis image; colormap(gray(256));
subplot(2,2,3:4); imagesc(mass_bg64); axis image; colormap(gray(256));

figure;
subplot(2,2,1); imagesc(norm_bg); axis image; colormap(gray(256));
subplot(2,2,2); imagesc(mass_bg); axis image; colormap(gray(256));
subplot(2,2,3:4); imagesc(mass_bg90); axis image; colormap(gray(256));
%%
figure; imagesc(norm_bg); axis image; colormap(gray(256)); caxis([min(norm_bg(:)) max(norm_bg(:))]);
figure; imagesc(mass_bg); axis image; colormap(gray(256));caxis([min(norm_bg(:)) max(norm_bg(:))]);
figure; imagesc(mass_bg64); axis image; colormap(gray(256));caxis([min(norm_bg(:)) max(norm_bg(:))]);
figure; imagesc(mass_bg90); axis image; colormap(gray(256));caxis([min(norm_bg(:)) max(norm_bg(:))]);
%%
write_im_from_colormap(norm_bg,'C:\isbe\thesis\figures\8\rot_norm_bg.bmp', colormap(gray(256)),[min(norm_bg(:)) max(norm_bg(:))]);
write_im_from_colormap(mass_bg,'C:\isbe\thesis\figures\8\rot_mass_bg.bmp', colormap(gray(256)),[min(norm_bg(:)) max(norm_bg(:))]);
write_im_from_colormap(mass_bg64,'C:\isbe\thesis\figures\8\rot_mass_bg64.bmp', colormap(gray(256)),[min(norm_bg(:)) max(norm_bg(:))]);
write_im_from_colormap(mass_bg90,'C:\isbe\thesis\figures\8\rot_mass_bg90.bmp', colormap(gray(256)),[min(norm_bg(:)) max(norm_bg(:))]);
%%

[bg64_orientation_map] = ...
        mb_dt_orientation_map('Image', mass_bg64);
[tout64 rout64 target_orientation_hist] = ...
        weighted_complex_rose(bg64_orientation_map(:), 120);
    
[bg90_orientation_map] = ...
        mb_dt_orientation_map('Image', mass_bg90);
[tout90 rout90 target_orientation_hist] = ...
        weighted_complex_rose(bg90_orientation_map(:), 120);
    
figure;
subplot(1,2,1); polar(tout64, rout64); title('BG 64 orientation');
subplot(1,2,2); polar(tout90, rout90); title('BG 90 orientation');
%%
max_radius = 140;
[x y] = meshgrid(1:size(norm_bg,2),(1:size(norm_bg,1))'); 
    target_mask = (x-384).^2 + (y-384).^2 < max_radius^2;
    
[norm_orientation_map] = ...
        mb_dt_orientation_map('Image', norm_bg);
[toutn routn target_orientation_hist] = ...
        weighted_complex_rose(norm_orientation_map(:), 120);
    
[mass_orientation_map] = ...
        mb_dt_orientation_map('Image', mass_bg);
[toutm routm target_orientation_hist] = ...
        weighted_complex_rose(mass_orientation_map(:), 120);
    
figure;
subplot(1,2,1); polar(toutn, routn); title('Target orientation');
subplot(1,2,2); polar(toutm, routm); title('Template orientation');
%%
f1 = figure(...
    'windowstyle', 'normal',...
    'Units', 'pixels',...
    'position', [100 100 450 450],...
    'PaperPositionMode','auto');
polar(toutn, routn);
print('-dtiff', '-noui', '-painters', f1, '-r300', 'C:\isbe\thesis\figures\8\rot_norm_bg_hist.tif');
f2 = figure(...
    'windowstyle', 'normal',...
    'Units', 'pixels',...
    'position', [100 100 450 450],...
    'PaperPositionMode','auto');
polar(toutm, routm);
print('-dtiff', '-noui', '-painters', f2, '-r300', 'C:\isbe\thesis\figures\8\rot_mass_bg_hist.tif');
%%
for r = 0:30:150
    [mass_rot] = ...
        mb_synthesise_mass(target_region, [512 512], mass_model,...
        mass_dir, model_id.mass_files, mm, r);
    figure; imagesc(mass_rot); axis image; colormap(gray(256));
end
%%
norm_list = dir('C:\isbe\dev\background\images\normal1024\*.bmp');
load C:\isbe\dev\mass_model\models\model_w500_50K.mat
mass_dir = 'C:\isbe\dev\masses1024x1024';
mkdir C:\isbe\thesis\figures\8\bgs
mass_rand = randsample(1:101, 10);
for nn = [48 78 81]
    for mmm = 1:10
        mm = mass_rand(mmm);
        try
            target_region = double(imread(['C:\isbe\dev\background\images\normal1024\', norm_list(nn).name]));
            [mass] = ...
                mb_synthesise_mass(target_region, [512 512], mass_model,...
                mass_dir, model_id.mass_files, mm, []);
            mass_name = ['C:\isbe\thesis\figures\8\bgs\norm', zerostr(nn,3), 'mass', zerostr(mm,3), '.bmp'];
            write_im_from_colormap(mass, mass_name, colormap(gray(256)));
        catch
            display('Skipped');
        end
    end
%     figure;
%     subplot(1,2,1); imagesc(target_region); axis image; colormap(gray(256));
%     subplot(1,2,2); imagesc(mass.mass_ROI); axis image; colormap(gray(256));
    
end
%%
load C:\isbe\dev\mass_model\models\model_w500_50K.mat
mass_dir = 'C:\isbe\dev\masses1024x1024';
norm = double(imread('C:\isbe\dev\background\images\normal1024\o04_023LCC_1024_3549_938.bmp'));
for mm = 1%:40
    try
        [mass] = ...
            mb_synthesise_mass(norm, [512 512], mass_model,...
            mass_dir, model_id.mass_files, mm, []);
        mass_name = ['C:\isbe\thesis\figures\8\norm_43_mass', zerostr(mm,3), '.bmp'];
        write_im_from_colormap(mass, mass_name, colormap(gray(256)));
    catch
        display('Skipped');
    end
end
%%
load('C:\isbe\dev\background\syn\dual_tree\normal512_normal512_k10_w1_5_w2_1_a_012.mat')
syn_name = ['C:\isbe\thesis\figures\8\syn_norm012.bmp'];
write_im_from_colormap(synthesised_image, syn_name, colormap(gray(256)));
load('C:\isbe\dev\background\syn\dual_tree\normal512_normal512_k10_w1_5_w2_1_a_013.mat')
syn_name = ['C:\isbe\thesis\figures\8\syn_norm013.bmp'];
write_im_from_colormap(synthesised_image, syn_name, colormap(gray(256)));
load('C:\isbe\dev\background\syn\dual_tree\normal512_normal512_k10_w1_5_w2_1_a_021.mat')
syn_name = ['C:\isbe\thesis\figures\8\syn_norm021.bmp'];
write_im_from_colormap(synthesised_image, syn_name, colormap(gray(256)));
load('C:\isbe\dev\background\syn\dual_tree\normal512_normal512_k10_w1_5_w2_1_a_043.mat')
syn_name = ['C:\isbe\thesis\figures\8\syn_norm043.bmp'];
write_im_from_colormap(synthesised_image, syn_name, colormap(gray(256)));
%%
mass_list = dir('C:\isbe\dev\masses1024x1024\*.mat');
for ii = [28 46 47]
    mass = u_load(['C:\isbe\dev\masses1024x1024\', mass_list(ii).name]);
    figure; imagesc(mass.mass_ROI); axis image; colormap(gray(256));
    clear mass;
end
%%
mass_list = dir('C:\isbe\dev\masses1024x1024\*.mat');
mass = u_load(['C:\isbe\dev\masses1024x1024\', mass_list(46).name]);
write_im_from_colormap(mass.mass_ROI, 'C:\isbe\thesis\figures\8\s_o_w.bmp',colormap(gray(256)));
write_im_from_colormap(mass.mass_ROI(350:650,500:800), 'C:\isbe\thesis\figures\8\s_o_w_zoom.bmp',colormap(gray(256)));