%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%% Make orientation images:
% mkdir C:\isbe\asymmetry_project\data\synthetic_lines\ori_syn512\labels
% mkdir C:\isbe\asymmetry_project\data\synthetic_lines\ori_real512\labels
mkdir C:\isbe\asymmetry_project\data\synthetic_lines\ori_real512_90\labels

for ii = [1 31 61]
    
%     bg_syn = u_load(['C:\isbe\asymmetry_project\data\synthetic_backgrounds\smooth512\test\bg' zerostr(ii,5) '.mat']);
    bg_real = u_load(['C:\isbe\asymmetry_project\data\synthetic_backgrounds\real512\test\bg' zerostr(ii,5) '.mat']);
    
    for jj = 1:6:180
        
        [bar_image, label, label_centre, label_orientation] = create_sin_bar(...
            3, 4, jj, 512, 512, 0.5, 256, 256);
        
%         bar_syn = bar_image + bg_syn;
%         bar_real = bar_image + bg_real;
        bar_real_90 = bar_image + rot90(bg_real);
        
%         save(['C:\isbe\asymmetry_project\data\synthetic_lines\ori_syn512\image'...
%             zerostr(jj,3) '_' zerostr(ii,3) '.mat'], 'bar_syn');
%         save(['C:\isbe\asymmetry_project\data\synthetic_lines\ori_real512\image'...
%             zerostr(jj,3) '_' zerostr(ii,3) '.mat'], 'bar_real');
        save(['C:\isbe\asymmetry_project\data\synthetic_lines\ori_real512_90\image'...
            zerostr(jj,3) '_' zerostr(ii,3) '.mat'], 'bar_real_90');
        
%         save(['C:\isbe\asymmetry_project\data\synthetic_lines\ori_syn512\labels\label'...
%             zerostr(jj,3) '_' zerostr(ii,3) '.mat'],...
%             'label', 'label_centre', 'label_orientation');
%         save(['C:\isbe\asymmetry_project\data\synthetic_lines\ori_real512\labels\label'...
%             zerostr(jj,3) '_' zerostr(ii,3) '.mat'],...
%             'label', 'label_centre', 'label_orientation');
        save(['C:\isbe\asymmetry_project\data\synthetic_lines\ori_real512_90\labels\label'...
            zerostr(jj,3) '_' zerostr(ii,3) '.mat'],...
            'label', 'label_centre', 'label_orientation');
    end
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%% Make width images:
mkdir C:\isbe\asymmetry_project\data\synthetic_lines\width_syn512\labels
mkdir C:\isbe\asymmetry_project\data\synthetic_lines\width_real512\labels

for ii = [1 31 61]
    
    bg_syn = u_load(['C:\isbe\asymmetry_project\data\synthetic_backgrounds\smooth512\test\bg' zerostr(ii,5) '.mat']);
    bg_real = u_load(['C:\isbe\asymmetry_project\data\synthetic_backgrounds\real512\test\bg' zerostr(ii,5) '.mat']);
    
    for jj = 1:15
        
        [bar_image, label, label_centre, label_orientation] = create_sin_bar(...
            jj, 4, 23, 512, 512, 0.5, 256, 256);
        
        bar_syn = bar_image + bg_syn;
        bar_real = bar_image + bg_real;
        
        save(['C:\isbe\asymmetry_project\data\synthetic_lines\width_syn512\image'...
            zerostr(jj,2) '_' zerostr(ii,3) '.mat'], 'bar_syn');
        save(['C:\isbe\asymmetry_project\data\synthetic_lines\width_real512\image'...
            zerostr(jj,2) '_' zerostr(ii,3) '.mat'], 'bar_real');
        
        save(['C:\isbe\asymmetry_project\data\synthetic_lines\width_syn512\labels\label'...
            zerostr(jj,2) '_' zerostr(ii,3) '.mat'],...
            'label', 'label_centre', 'label_orientation');
        save(['C:\isbe\asymmetry_project\data\synthetic_lines\width_real512\labels\label'...
            zerostr(jj,2) '_' zerostr(ii,3) '.mat'],...
            'label', 'label_centre', 'label_orientation');
    end
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Make contrast images:
mkdir C:\isbe\asymmetry_project\data\synthetic_lines\con_syn512\labels
mkdir C:\isbe\asymmetry_project\data\synthetic_lines\con_real512\labels

for ii = [1 31 61]
    
    bg_syn = u_load(['C:\isbe\asymmetry_project\data\synthetic_backgrounds\smooth512\test\bg' zerostr(ii,5) '.mat']);
    bg_real = u_load(['C:\isbe\asymmetry_project\data\synthetic_backgrounds\real512\test\bg' zerostr(ii,5) '.mat']);
    
    for jj = 2:16
        
        [bar_image, label, label_centre, label_orientation] = create_sin_bar(...
            3, jj, 23, 512, 512, 0.5, 256, 256);
        
        bar_syn = bar_image + bg_syn;
        bar_real = bar_image + bg_real;
        
        save(['C:\isbe\asymmetry_project\data\synthetic_lines\con_syn512\image'...
            zerostr(jj,2) '_' zerostr(ii,3) '.mat'], 'bar_syn');
        save(['C:\isbe\asymmetry_project\data\synthetic_lines\con_real512\image'...
            zerostr(jj,2) '_' zerostr(ii,3) '.mat'], 'bar_real');
        
        save(['C:\isbe\asymmetry_project\data\synthetic_lines\con_syn512\labels\label'...
            zerostr(jj,2) '_' zerostr(ii,3) '.mat'],...
            'label', 'label_centre', 'label_orientation');
        save(['C:\isbe\asymmetry_project\data\synthetic_lines\con_real512\labels\label'...
            zerostr(jj,2) '_' zerostr(ii,3) '.mat'],...
            'label', 'label_centre', 'label_orientation');
    end
end        
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Make squash images:
mkdir C:\isbe\asymmetry_project\data\synthetic_lines\squash_syn512\labels
mkdir C:\isbe\asymmetry_project\data\synthetic_lines\squash_real512\labels

for ii = [1 31 61]
    
    bg_syn = u_load(['C:\isbe\asymmetry_project\data\synthetic_backgrounds\smooth512\test\bg' zerostr(ii,5) '.mat']);
    bg_real = u_load(['C:\isbe\asymmetry_project\data\synthetic_backgrounds\real512\test\bg' zerostr(ii,5) '.mat']);
    
    for jj = 0:10
        
        [bar_image, label, label_centre, label_orientation] = create_sin_bar(...
            3, 4, 23, 512, 512, jj/10, 256, 256);
        
        bar_syn = bar_image + bg_syn;
        bar_real = bar_image + bg_real;
        
        save(['C:\isbe\asymmetry_project\data\synthetic_lines\squash_syn512\image'...
            zerostr(jj,2) '_' zerostr(ii,3) '.mat'], 'bar_syn');
        save(['C:\isbe\asymmetry_project\data\synthetic_lines\squash_real512\image'...
            zerostr(jj,2) '_' zerostr(ii,3) '.mat'], 'bar_real');
        
        save(['C:\isbe\asymmetry_project\data\synthetic_lines\squash_syn512\labels\label'...
            zerostr(jj,2) '_' zerostr(ii,3) '.mat'],...
            'label', 'label_centre', 'label_orientation');
        save(['C:\isbe\asymmetry_project\data\synthetic_lines\squash_real512\labels\label'...
            zerostr(jj,2) '_' zerostr(ii,3) '.mat'],...
            'label', 'label_centre', 'label_orientation');
    end
end 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%% Compute orientation maps for g2 method (RF computed on hydra)
%--------------------------------------------------------------------------
% Notes: substitute 'real' for 'syn'
%% Orientations

image_type = 'real';
mkdir(['C:\isbe\asymmetry_project\data\synthetic_lines\ori_' image_type '512_90\results\g2_ori']);

for jj = 1:6:180
    for ii = [1 31 61]
        
        test_image = load_uint8(['C:\isbe\asymmetry_project\data\synthetic_lines\ori_' image_type '512_90\image'...
            zerostr(jj,3) '_' zerostr(ii,3) '.mat']);
        
        [orientation_map] = karssemeijer_line_detection(...
            test_image,...
            'line_scales', [1 2 4 8],...
            'grad_scale', 10,...
            'grad_ori_thresh', pi/6,...
            'grad_strength_thresh', 25,...
            'line_strength_thresh', 0,...
            'binary_map', 1);
        
        save_uint8(['C:\isbe\asymmetry_project\data\synthetic_lines\ori_' image_type '512\results\g2_ori\image'...
            zerostr(jj,3) '_' zerostr(ii,3) '_class.mat'],  orientation_map);
        
    end
end
%%
% Widths
mkdir(['C:\isbe\asymmetry_project\data\synthetic_lines\width_' image_type '512\results\g2_ori']);

for jj = 1:15
    for ii = [1 31 61]
        
        test_image = load_uint8(['C:\isbe\asymmetry_project\data\synthetic_lines\width_' image_type '512\image'...
            zerostr(jj,2) '_' zerostr(ii,3) '.mat']);
        
        [orientation_map] = karssemeijer_line_detection(...
            test_image,...
            'line_scales', [1 2 4 8],...
            'grad_scale', 10,...
            'grad_ori_thresh', pi/6,...
            'grad_strength_thresh', 25,...
            'line_strength_thresh', 0,...
            'binary_map', 1);
        
        save_uint8(['C:\isbe\asymmetry_project\data\synthetic_lines\width_' image_type '512\results\g2_ori\image'...
            zerostr(jj,2) '_' zerostr(ii,3) 'ori.mat'],  orientation_map);
        
    end
end
% Contrast
mkdir(['C:\isbe\asymmetry_project\data\synthetic_lines\con_' image_type '512\results\g2_ori']);
for jj = 2:16
    for ii = [1 31 61]
        
        test_image = load_uint8(['C:\isbe\asymmetry_project\data\synthetic_lines\con_' image_type '512\image'...
            zerostr(jj,2) '_' zerostr(ii,3) '.mat']);
        
        [orientation_map] = karssemeijer_line_detection(...
            test_image,...
            'line_scales', [1 2 4 8],...
            'grad_scale', 10,...
            'grad_ori_thresh', pi/6,...
            'grad_strength_thresh', 25,...
            'line_strength_thresh', 0,...
            'binary_map', 1);
        
        save_uint8(['C:\isbe\asymmetry_project\data\synthetic_lines\con_' image_type '512\results\g2_ori\image'...
            zerostr(jj,2) '_' zerostr(ii,3) 'ori.mat'],  orientation_map);
        
    end
end

% Squash
mkdir(['C:\isbe\asymmetry_project\data\synthetic_lines\squash_' image_type '512\results\g2_ori']);
for ii = [1 31 61]

    for jj = 0:10
        
        test_image = load_uint8(['C:\isbe\asymmetry_project\data\synthetic_lines\squash_' image_type '512\image'...
            zerostr(jj,2) '_' zerostr(ii,3) '.mat']);
        
        [orientation_map] = karssemeijer_line_detection(...
            test_image,...
            'line_scales', [1 2 4 8],...
            'grad_scale', 10,...
            'grad_ori_thresh', pi/6,...
            'grad_strength_thresh', 25,...
            'line_strength_thresh', 0,...
            'binary_map', 1);
        
        save_uint8(['C:\isbe\asymmetry_project\data\synthetic_lines\squash_' image_type '512\results\g2_ori\image'...
            zerostr(jj,2) '_' zerostr(ii,3) 'ori.mat'],  orientation_map);

    end
end
%--------------------------------------------------------------------------
copyfile('Z:\data\synthetic_lines\ori_real512\results\238470',...
    'C:\isbe\asymmetry_project\data\synthetic_lines\ori_real512\results\238470');

copyfile('Z:\data\synthetic_lines\con_real512\results\238470',...
    'C:\isbe\asymmetry_project\data\synthetic_lines\con_real512\results\238470');

copyfile('Z:\data\synthetic_lines\width_real512\results\238470',...
    'C:\isbe\asymmetry_project\data\synthetic_lines\width_real512\results\238470');

copyfile('Z:\data\synthetic_lines\squash_real512\results\238470',...
    'C:\isbe\asymmetry_project\data\synthetic_lines\squash_real512\results\238470');
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%% Produce error plots of the various image for the g2/RF methods
%--------------------------------------------------------------------------
%% Orientations
image_type = 'real';

xy = repmat(-256:255, 512, 1);
circle_mask = xy.^2 + xy'.^2 < 128^2;
figure; hold on;
mean_errors = [];

for jj = 1:6:180
    
    ori_errors = [];
    for ii = [1 31 61]
        
        ori_map_rf = load_uint8(['C:\isbe\asymmetry_project\data\synthetic_lines\ori_' image_type '512\results\238470\image'...
            zerostr(jj,3) '_' zerostr(ii,3) '_class.mat']);
        ori_map_g2 = load_uint8(['C:\isbe\asymmetry_project\data\synthetic_lines\ori_' image_type '512\results\g2_ori\image'...
            zerostr(jj,3) '_' zerostr(ii,3) 'ori.mat']);
        label = load(['C:\isbe\asymmetry_project\data\synthetic_lines\ori_' image_type '512\labels\label'...
            zerostr(jj,3) '_' zerostr(ii,3) '.mat']);

        label_centre = label.label_centre & circle_mask;
        
        ori_errors = [ori_errors; ...
            abs(mb_mod(ori_map_rf(label_centre) - label.label_orientation(label_centre), 180))...
            abs(mb_mod(ori_map_g2(label_centre) - label.label_orientation(label_centre), 180))];%#ok
        %plot(jj, mean(abs(ori_errors)), 'x');
        display(['Angle = ' num2str(jj)]);
        display(min(ori_errors));
        display(max(ori_errors));
        
    end
    boxplot(ori_errors, 'positions', jj+[-1 1], 'labels', {'',''}, 'widths', 2, 'whisker', 2, 'colors', 'kb');
    %set(h(7), 'visible', 'off');
    mean_errors(end+1,:) = mean(ori_errors); %#ok
end
plot(1:6:180, mean_errors(:,1), 'g', 'linewidth', 2);
plot(1:6:180, mean_errors(:,2), 'c', 'linewidth', 2);
axis([-4 180 0 90]);
title('Orientation errors for lines of varying orientation');
ylabel('Mean absolute orientation error (degrees)');
xlabel('Line orientation (degrees)');
set(gca, 'xtick', 1:6:180, 'xticklabel', 1:6:180);
%% Widths
figure; hold on;
mean_errors = zeros(15,2);
for jj = 1:15
    
    ori_errors = [];
    for ii = [1 31 61]
        
        ori_map_rf = load_uint8(['C:\isbe\asymmetry_project\data\synthetic_lines\width_' image_type '512\results\238470\image'...
            zerostr(jj,2) '_' zerostr(ii,3) '_class.mat']);
        ori_map_g2 = load_uint8(['C:\isbe\asymmetry_project\data\synthetic_lines\width_' image_type '512\results\g2_ori\image'...
            zerostr(jj,2) '_' zerostr(ii,3) 'ori.mat']);
        label = load(['C:\isbe\asymmetry_project\data\synthetic_lines\width_' image_type '512\labels\label'...
            zerostr(jj,2) '_' zerostr(ii,3) '.mat']);

        ori_errors = [ori_errors; ...
            abs(mb_mod(ori_map_rf(label.label_centre) - label.label_orientation(label.label_centre), 180))...
            abs(mb_mod(ori_map_g2(label.label_centre) - label.label_orientation(label.label_centre), 180))];%#ok
        %plot(jj, mean(abs(ori_errors)), 'x');
        
    end
    boxplot(ori_errors, 'positions', jj+[-.2 .2], 'labels', {'',''}, 'widths', 0.4, 'whisker', 2, 'colors', 'kb');
    mean_errors(jj,:) = mean(ori_errors);

end
plot(1:15, mean_errors(:,1), 'g', 'linewidth', 2);
plot(1:15, mean_errors(:,2), 'c', 'linewidth', 2);
axis([-.5 16 0 90]);
title('Orientation errors for lines of varying width');
ylabel('Mean absolute orientation error');
xlabel('Line half-width (pixels)');
set(gca, 'xtick', 1:15, 'xticklabel', 1:15);
%% Contrast
figure; hold on;
mean_errors = zeros(15,2);
for jj = 2:16
    ori_errors = [];
    for ii = [1 31 61]
        
        ori_map_rf = load_uint8(['C:\isbe\asymmetry_project\data\synthetic_lines\con_' image_type '512\results\238470\image'...
            zerostr(jj,2) '_' zerostr(ii,3) '_class.mat']);
        ori_map_g2 = load_uint8(['C:\isbe\asymmetry_project\data\synthetic_lines\con_' image_type '512\results\g2_ori\image'...
            zerostr(jj,2) '_' zerostr(ii,3) 'ori.mat']);
        label = load(['C:\isbe\asymmetry_project\data\synthetic_lines\con_' image_type '512\labels\label'...
            zerostr(jj,2) '_' zerostr(ii,3) '.mat']);

        ori_errors = [ori_errors; ...
            abs(mb_mod(ori_map_rf(label.label_centre) - label.label_orientation(label.label_centre), 180))...
            abs(mb_mod(ori_map_g2(label.label_centre) - label.label_orientation(label.label_centre), 180))];%#ok
        %plot(jj, mean(abs(ori_errors)), 'x');
        
    end
    boxplot(ori_errors, 'positions', jj+[-.2 .2], 'labels', {'',''}, 'widths', 0.4, 'whisker', 2, 'colors', 'kb');
    mean_errors(jj-1,:) = mean(ori_errors);
end
plot(2:16, mean_errors(:,1), 'g', 'linewidth', 2);
plot(2:16, mean_errors(:,2), 'c', 'linewidth', 2);
axis([0 17 0 90]);
title('Orientation errors for lines of varying line contrast');
ylabel('Mean absolute orientation error');
xlabel('Contrast at line centre (8-bit grey-level)');
set(gca, 'xtick', 2:16, 'xticklabel', 2:16);
%% Squash
figure; hold on;
mean_errors = zeros(11,2);
for jj = 0:10
    ori_errors = [];
    
    for ii = [1 31 61]
        ori_map_rf = load_uint8(['C:\isbe\asymmetry_project\data\synthetic_lines\squash_' image_type '512\results\238470\image'...
            zerostr(jj,2) '_' zerostr(ii,3) '_class.mat']);
        ori_map_g2 = load_uint8(['C:\isbe\asymmetry_project\data\synthetic_lines\squash_' image_type '512\results\g2_ori\image'...
            zerostr(jj,2) '_' zerostr(ii,3) 'ori.mat']);
        label = load(['C:\isbe\asymmetry_project\data\synthetic_lines\squash_' image_type '512\labels\label'...
            zerostr(jj,2) '_' zerostr(ii,3) '.mat']);

        ori_errors = [ori_errors; ...
            abs(mb_mod(ori_map_rf(label.label_centre) - label.label_orientation(label.label_centre), 180))...
            abs(mb_mod(ori_map_g2(label.label_centre) - label.label_orientation(label.label_centre), 180))];%#ok
        %plot(jj, mean(abs(ori_errors)), 'x');
        
    end
    boxplot(ori_errors, 'positions', jj+[-.2 .2], 'labels', {'',''}, 'widths', 0.4, 'whisker', 2, 'colors', 'kb');
    mean_errors(jj+1,:) = mean(ori_errors);
end
plot(0:10, mean_errors(:,1), 'g', 'linewidth', 2);
plot(0:10, mean_errors(:,2), 'c', 'linewidth', 2);
axis([-1 11 0 90]);
title('Orientation errors for lines of varying squashiness');
ylabel('Mean absolute orientation error');
xlabel('Squashiness (0 - 1)');
set(gca, 'xtick', 0:10, 'xticklabel', 0:10);
%%
angles = [];
widths = [];
contrasts = [];
squashes = [];
for rf = 1:20
    for tree = 1:10
        par = u_load(['Z:\data\line_orientation_rfs\233902\line_parameters\' zerostr(rf,2) '\parameters' zerostr(tree,3)]);
        for ii = 1:length(par)
            angles(end+1,1) = mod(par(ii).orientation,180);
            widths(end+1,1) = par(ii).width;
            contrasts(end+1,1) = par(ii).contrast;
            squashes(end+1,1) = par(ii).squash;
        end
    end
end
%%
angles = [];
widths = [];
contrasts = [];
squashes = [];
for rf = 1:20
    for tree = 1:10
        par = u_load(['Z:\data\line_orientation_rfs\233902\line_parameters\' zerostr(rf,2) '\parameters' zerostr(tree,3)]);
        for ii = 1:length(par)
            if ismember(par(ii).bg_idx, [1 31 61])
                angles(end+1,1) = mod(par(ii).orientation,180);
                widths(end+1,1) = par(ii).width;
                contrasts(end+1,1) = par(ii).contrast;
                squashes(end+1,1) = par(ii).squash;
            end
        end
    end
end
%%
figure; hist(angles, 180); title('Histogram of training data orientations');
figure; hist(widths/2, 100); title('Histogram of training data widths');
figure; hist(contrasts, 100); title('Histogram of training data contrasts');
figure; hist(squashes, 100); title('Histogram of training data squashes');
%%
mkdir C:\isbe\asymmetry_project\data\synthetic_lines\ori_real512_90\
im_list = dir('C:\isbe\asymmetry_project\data\synthetic_lines\ori_real512\*.mat');
for ii = 1:length(im_list)
    test_im = u_load(['C:\isbe\asymmetry_project\data\synthetic_lines\ori_real512\' im_list(ii).name]);
    test_im = rot90(test_im);
    save(['C:\isbe\asymmetry_project\data\synthetic_lines\ori_real512_90\' im_list(ii).name], 'test_im');
end
%%
image_type = 'real';

figure; hold on;
mean_errors = [];
ssq_errors = [];

for jj = 1:6:180
    
    ori_errors = [];
    for ii = [1 51 101]
        
        ori_map_rf = load_uint8(['C:\isbe\asymmetry_project\data\synthetic_lines\old_lines\ori_' image_type '512_90\results\233902\image'...
            zerostr(jj,3) '_' zerostr(ii,3) '_class.mat']);
        ori_map_g2 = load_uint8(['C:\isbe\asymmetry_project\data\synthetic_lines\old_lines\ori_' image_type '512\results\g2_ori\image'...
           zerostr(jj,3) '_' zerostr(ii,3) 'ori.mat']);
%         ori_map_g2 = load_uint8(['C:\isbe\asymmetry_project\data\synthetic_lines\ori_' image_type '512\results\233902\image'...
%             zerostr(jj,3) '_' zerostr(ii,3) '_class.mat']);
        label = load(['C:\isbe\asymmetry_project\data\synthetic_lines\old_lines\ori_' image_type '512\labels\label'...
            zerostr(jj,3) '_' zerostr(ii,3) '.mat']);

        lab_ori = rot90(label.label_orientation)+90;
        lab_centre = rot90(label.label_centre);
        ori_errors = [ori_errors; ...
            abs(mb_mod(ori_map_rf(lab_centre) - lab_ori(lab_centre), 180))...
            abs(mb_mod(ori_map_g2(label.label_centre) - label.label_orientation(label.label_centre), 180))];%#ok
        %plot(jj, mean(abs(ori_errors)), 'x');
        
    end
    boxplot(ori_errors, 'positions', jj+[-1 1], 'labels', {'',''}, 'widths', 2, 'whisker', 2, 'colors', 'kb');
    %set(h(7), 'visible', 'off');
    mean_errors(end+1,:) = mean(ori_errors); %#ok
end
plot(1:6:180, mean_errors(:,1), 'g', 'linewidth', 2);
plot(1:6:180, mean_errors(:,2), 'c', 'linewidth', 2);
axis([-1 180 0 90]);
title('Orientation errors for lines of varying orientation');
ylabel('Mean absolute orientation error (degrees)');
xlabel('Line orientation (degrees)');
set(gca, 'xtick', 1:6:180, 'xticklabel', 1:6:180);
%%
for jj = [1 7 13 85 91 97 169 175]
    for ii = [1 31 61]
        
        %figure; imagesc(u_load(['C:\isbe\asymmetry_project\data\synthetic_lines\ori_' image_type '512_90\image'...
        %    zerostr(jj,3) '_' zerostr(ii,3) '.mat'])); axis image; colormap(gray(256));
        figure; imagesc(u_load(['C:\isbe\asymmetry_project\data\synthetic_lines\ori_' image_type '512_90\results\233902\image'...
            zerostr(jj,3) '_' zerostr(ii,3) '_class.mat'])); axis image; colormap(hsv(180));
    end
end
%%
image_type = 'real';

for ii = [1 31 61]
    mean_errors = [];
    figure; hold on;
    for jj = 1:6:180
    
 
        
        ori_map_rf = load_uint8(['C:\isbe\asymmetry_project\data\synthetic_lines\ori_' image_type '512\results\233902\image'...
            zerostr(jj,3) '_' zerostr(ii,3) '_class.mat']);
        ori_map_g2 = load_uint8(['C:\isbe\asymmetry_project\data\synthetic_lines\ori_' image_type '512\results\g2_ori\image'...
            zerostr(jj,3) '_' zerostr(ii,3) 'ori.mat']);
        label = load(['C:\isbe\asymmetry_project\data\synthetic_lines\ori_' image_type '512\labels\label'...
            zerostr(jj,3) '_' zerostr(ii,3) '.mat']);

        %lab_ori = rot90(label.label_orientation)-90;
        lab_ori = label.label_orientation;
        
        ori_errors = [...
            abs(mb_mod(ori_map_rf(label.label_centre) - lab_ori(label.label_centre), 180))...
            abs(mb_mod(ori_map_g2(label.label_centre) - label.label_orientation(label.label_centre), 180))];
        %plot(jj, mean(abs(ori_errors)), 'x');
        boxplot(ori_errors, 'positions', jj+[-1 1], 'labels', {'',''}, 'widths', 2, 'whisker', 2, 'colors', 'kb');
        mean_errors(end+1,:) = mean(ori_errors); %#ok
    end
    plot(1:6:180, mean_errors(:,1), 'g', 'linewidth', 2);
    plot(1:6:180, mean_errors(:,2), 'c', 'linewidth', 2);
    axis([-1 180 0 90]);
    title('Orientation errors for lines of varying orientation');
    ylabel('Mean absolute orientation error (degrees)');
    xlabel('Line orientation (degrees)');
    set(gca, 'xtick', 1:6:180, 'xticklabel', 1:6:180);   
end
%%
% Just checking like, but are we sure the 180 degree wrapping in the
% regression prediction is working...

load C:\isbe\asymmetry_project\data\misc\y_trees y_trees
for jj = -90:10:90
    y_fit = 0;
    y_tree = y_trees(121,:) + jj;
    for ii = 1:200

        y_fit = mod(y_fit - mb_mod(y_fit - y_tree(ii),180)/ii,180);
    end
    display([num2str(mod(y_fit - jj,180)) '   ' num2str(y_fit)]);
end
%%
for jj = 1:10
    y_fit = 0;
    y_tree = y_trees(121,:);
    y_tree = y_tree(randperm(200));
    for ii = 1:200

        y_fit = mod(y_fit - mb_mod(y_fit - y_tree(ii),180)/ii,180);
    end
    display([num2str(y_fit) '   ' num2str(mean(y_tree))]);
end
%%
clc;
n = 6;
pp = perms(1:n);
for kk = 0:(200-n)
    y_fits = zeros(factorial(n),1);
    for jj = 1:factorial(n)
        y_fit = 0;
        y_tree = y_trees(121,kk+(1:n));
        y_tree = y_tree(pp(jj,:));
        for ii = 1:n

            %y_fit = mod(y_fit - mb_mod(y_fit - y_tree(ii),180)/ii,180);
            %y_fit = mod(y_fit - (y_fit - y_tree(ii))/ii,180);
            y_fit = y_fit - mb_mod(y_fit - y_tree(ii),180)/ii;
            
            %display(num2str(y_fit));
        end
        y_fits(jj) = y_fit;
        %display('******');
    end
    if any(diff(y_fits)>1e-6)
        display(kk);
    end
end
%%
clc;
n = 5;
pp = perms(1:n);
kk = 1;
y_fits = zeros(factorial(n),1);
for jj = 2:3%1:factorial(n)
    y_fit = 0;
    y_tree = y_trees(121,kk+(1:n));
    y_tree = y_tree(pp(jj,:));
    
    figure;
    for ii = 1:n
        subplot(2,3,ii); axis([-1 1 -1 1]); hold on;
        plot([0 cosd(y_fit)], [0 sind(y_fit)], 'b');
        plot([0 cosd(y_tree(ii))], [0 sind(y_tree(ii))], 'b');
        
        y_fit = mod(y_fit - mb_mod(y_fit - y_tree(ii),180)/ii,180);
        display(num2str(y_fit));
        
        
    end
    y_fits(jj) = y_fit;
    display('******');
end
%%
load C:\isbe\asymmetry_project\data\misc\y_trees y_trees

[N ntrees] = size(y_trees);
y_fit = zeros(N,1);

for ii = 1:ntrees
    y_fit = mod(y_fit - mb_mod(y_fit - y_trees(:,ii),180)/ii,180);
end
%%
ori_estimates = pi*y_trees/90; %[0, 2pi]
x_mean = mean(cos(ori_estimates),2);
y_mean = mean(sin(ori_estimates),2);
y_fit2 = mod(90*atan2(y_mean, x_mean)/pi,180);
y_r = sqrt(x_mean.^2 + y_mean.^2);
%%
image_type = 'real';

xy = repmat(-256:255, 512, 1);
circle_mask = xy.^2 + xy'.^2 < 128^2;
figure; hold on;
mean_errors = [];

for jj = 1:6:180
    
    ori_errors = [];
    for ii = [1 31 61]
        
        ori_map_rf1 = mod(180*angle(load_uint8(['C:\isbe\asymmetry_project\data\synthetic_lines\ori_' image_type '512\results\213208\image'...
            zerostr(jj,3) '_' zerostr(ii,3) '_class.mat']))/pi,180);
        ori_map_rf2 = mod(180*angle(load_uint8(['C:\isbe\asymmetry_project\data\synthetic_lines\ori_' image_type '512\results\238470\image'...
            zerostr(jj,3) '_' zerostr(ii,3) '_class.mat']))/pi,180);
        ori_map_g2 = load_uint8(['C:\isbe\asymmetry_project\data\synthetic_lines\ori_' image_type '512\results\g2_ori\image'...
           zerostr(jj,3) '_' zerostr(ii,3) '_class.mat']);
        label = load(['C:\isbe\asymmetry_project\data\synthetic_lines\ori_' image_type '512\labels\label'...
            zerostr(jj,3) '_' zerostr(ii,3) '.mat']);

        label_centre = label.label_centre & circle_mask;
        
        ori_errors = [ori_errors; ...
            abs(mb_mod(ori_map_rf1(label_centre) - label.label_orientation(label_centre), 180))...
            abs(mb_mod(ori_map_rf2(label_centre) - label.label_orientation(label_centre), 180))...
            abs(mb_mod(ori_map_g2(label_centre) - label.label_orientation(label_centre), 180))];%#ok
        %plot(jj, mean(abs(ori_errors)), 'x');
        display(['Angle = ' num2str(jj)]);
        display(min(ori_errors));
        display(max(ori_errors));
        
    end
    %boxplot(ori_errors, 'positions', jj+[-1 1], 'labels', {'',''}, 'widths', 2, 'whisker', 2, 'colors', 'kb');
    %set(h(7), 'visible', 'off');
    mean_errors(end+1,:) = mean(ori_errors); %#ok
end
plot(1:6:180, mean_errors(:,1), 'g', 'linewidth', 2);
plot(1:6:180, mean_errors(:,2), 'm', 'linewidth', 2);
plot(1:6:180, mean_errors(:,3), 'c', 'linewidth', 2);

legend({'RF trained on synthetic BGs', 'RF trained on real BGs', 'Gaussian derivatives'});

plot([1 180], [mean(mean_errors(:,1)) mean(mean_errors(:,1))], 'g--', 'linewidth', 2);
plot([1 180], [mean(mean_errors(:,2)) mean(mean_errors(:,2))], 'm--', 'linewidth', 2);
plot([1 180], [mean(mean_errors(:,3)) mean(mean_errors(:,3))], 'c--', 'linewidth', 2);

axis([-4 180 0 60]);
title('Orientation errors for lines of varying orientation');
ylabel('Mean absolute orientation error (degrees)');
xlabel('Line orientation (degrees)');
set(gca, 'xtick', 1:6:180, 'xticklabel', 1:6:180);
%%
image_type = 'real';

xy = repmat(-256:255, 512, 1);
circle_mask = xy.^2 + xy'.^2 < 128^2;
figure; hold on;
mean_errors = [];

for jj = 1:6:180
    
    ori_errors = [];
    for ii = [1 31 61]
        
        ori_map_rf1 = mod(180*angle(load_uint8(['C:\isbe\asymmetry_project\data\synthetic_lines\ori_' image_type '512_90\results\213208\image'...
            zerostr(jj,3) '_' zerostr(ii,3) '_class.mat']))/pi,180);
        ori_map_rf2 = mod(180*angle(load_uint8(['C:\isbe\asymmetry_project\data\synthetic_lines\ori_' image_type '512_90\results\238470\image'...
            zerostr(jj,3) '_' zerostr(ii,3) '_class.mat']))/pi,180);
        ori_map_g2 = load_uint8(['C:\isbe\asymmetry_project\data\synthetic_lines\ori_' image_type '512_90\results\g2_ori\image'...
           zerostr(jj,3) '_' zerostr(ii,3) '_class.mat']);
        label = load(['C:\isbe\asymmetry_project\data\synthetic_lines\ori_' image_type '512_90\labels\label'...
            zerostr(jj,3) '_' zerostr(ii,3) '.mat']);

        label_centre = label.label_centre & circle_mask;
        
        ori_errors = [ori_errors; ...
            abs(mb_mod(ori_map_rf1(label_centre) - label.label_orientation(label_centre), 180))...
            abs(mb_mod(ori_map_rf2(label_centre) - label.label_orientation(label_centre), 180))...
            abs(mb_mod(ori_map_g2(label_centre) - label.label_orientation(label_centre), 180))];%#ok
        %plot(jj, mean(abs(ori_errors)), 'x');
        display(['Angle = ' num2str(jj)]);
        display(min(ori_errors));
        display(max(ori_errors));
        
    end
    %boxplot(ori_errors, 'positions', jj+[-1 1], 'labels', {'',''}, 'widths', 2, 'whisker', 2, 'colors', 'kb');
    %set(h(7), 'visible', 'off');
    mean_errors(end+1,:) = mean(ori_errors); %#ok
end
plot(1:6:180, mean_errors(:,1), 'g', 'linewidth', 2);
plot(1:6:180, mean_errors(:,2), 'm', 'linewidth', 2);
plot(1:6:180, mean_errors(:,3), 'c', 'linewidth', 2);

legend({'RF trained on synthetic BGs', 'RF trained on real BGs', 'Gaussian derivatives'});

plot([1 180], [mean(mean_errors(:,1)) mean(mean_errors(:,1))], 'g--', 'linewidth', 2);
plot([1 180], [mean(mean_errors(:,2)) mean(mean_errors(:,2))], 'm--', 'linewidth', 2);
plot([1 180], [mean(mean_errors(:,3)) mean(mean_errors(:,3))], 'c--', 'linewidth', 2);

axis([-4 180 0 60]);
title('Orientation errors for lines of varying orientation');
ylabel('Mean absolute orientation error (degrees)');
xlabel('Line orientation (degrees)');
set(gca, 'xtick', 1:6:180, 'xticklabel', 1:6:180);
%%
image_type = {'real512', 'real512_90'};

xy = repmat(-256:255, 512, 1);
circle_mask = xy.^2 + xy'.^2 < 128^2;
figure; hold on;
mean_errors = [];

for jj = 1:6:180
    
    ori_errors = [];
    for kk = 1:2
        for ii = [1 31 61]

            ori_map_rf1 = mod(180*angle(load_uint8(['C:\isbe\asymmetry_project\data\synthetic_lines\ori_' image_type{kk} '\results\213208\image'...
                zerostr(jj,3) '_' zerostr(ii,3) '_class.mat']))/pi,180);
            ori_map_rf2 = mod(180*angle(load_uint8(['C:\isbe\asymmetry_project\data\synthetic_lines\ori_' image_type{kk} '\results\238470\image'...
                zerostr(jj,3) '_' zerostr(ii,3) '_class.mat']))/pi,180);
            ori_map_g2 = load_uint8(['C:\isbe\asymmetry_project\data\synthetic_lines\ori_' image_type{kk} '\results\g2_ori\image'...
               zerostr(jj,3) '_' zerostr(ii,3) '_class.mat']);
            label = load(['C:\isbe\asymmetry_project\data\synthetic_lines\ori_' image_type{kk} '\labels\label'...
                zerostr(jj,3) '_' zerostr(ii,3) '.mat']);

            label_centre = label.label_centre & circle_mask;

            ori_errors = [ori_errors; ...
                abs(mb_mod(ori_map_rf1(label_centre) - label.label_orientation(label_centre), 180))...
                abs(mb_mod(ori_map_rf2(label_centre) - label.label_orientation(label_centre), 180))...
                abs(mb_mod(ori_map_g2(label_centre) - label.label_orientation(label_centre), 180))];%#ok
            %plot(jj, mean(abs(ori_errors)), 'x');
            display(['Angle = ' num2str(jj)]);
            display(min(ori_errors));
            display(max(ori_errors));

        end
    end
    %boxplot(ori_errors, 'positions', jj+[-1 1], 'labels', {'',''}, 'widths', 2, 'whisker', 2, 'colors', 'kb');
    %set(h(7), 'visible', 'off');
    mean_errors(end+1,:) = mean(ori_errors); %#ok
end
plot(1:6:180, mean_errors(:,1), 'g', 'linewidth', 2);
plot(1:6:180, mean_errors(:,2), 'm', 'linewidth', 2);
plot(1:6:180, mean_errors(:,3), 'c', 'linewidth', 2);

legend({'RF trained on synthetic BGs', 'RF trained on real BGs', 'Gaussian derivatives'});

plot([1 180], [mean(mean_errors(:,1)) mean(mean_errors(:,1))], 'g--', 'linewidth', 2);
plot([1 180], [mean(mean_errors(:,2)) mean(mean_errors(:,2))], 'm--', 'linewidth', 2);
plot([1 180], [mean(mean_errors(:,3)) mean(mean_errors(:,3))], 'c--', 'linewidth', 2);

axis([-4 180 0 60]);
title('Orientation errors for lines of varying orientation');
ylabel('Mean absolute orientation error (degrees)');
xlabel('Line orientation (degrees)');
set(gca, 'xtick', 1:6:180, 'xticklabel', 1:6:180);
%%
image_type = 'real';

%xy = repmat(-256:255, 512, 1);
%circle_mask = xy.^2 + xy'.^2 < 128^2;
circle_mask = true(512);
f1 = figure;
f2 = figure;
f3 = figure;

count = 1;
for jj = 1:12:180
    
    ori_errors = [];
    for ii = [1 31 61]
        
        ori_map_rf1 = mod(180*angle(load_uint8(['C:\isbe\asymmetry_project\data\synthetic_lines\ori_' image_type '512_90\results\213208\image'...
            zerostr(jj,3) '_' zerostr(ii,3) '_class.mat']))/pi,180);
        ori_map_rf2 = mod(180*angle(load_uint8(['C:\isbe\asymmetry_project\data\synthetic_lines\ori_' image_type '512_90\results\238470\image'...
            zerostr(jj,3) '_' zerostr(ii,3) '_class.mat']))/pi,180);
        ori_map_g2 = mod(load_uint8(['C:\isbe\asymmetry_project\data\synthetic_lines\ori_' image_type '512_90\results\g2_ori\image'...
           zerostr(jj,3) '_' zerostr(ii,3) '_class.mat']),180);
        label = load(['C:\isbe\asymmetry_project\data\synthetic_lines\ori_' image_type '512_90\labels\label'...
            zerostr(jj,3) '_' zerostr(ii,3) '.mat']);

        label_centre = label.label_centre & circle_mask;
        
        ori_estimates = [ori_errors; ...
            ori_map_rf1(label_centre)...
            ori_map_rf2(label_centre)...
            ori_map_g2(label_centre)];%#ok
        
    end
    figure(f1); subplot(3,5,count); hist(ori_estimates(:,1), 0:180); hold on;
    plot([jj jj], get(gca, 'ylim'), 'r:');
    title([num2str(jj) '^{\circ}']); set(gca, 'xlim', [-1 181]);
    figure(f2); subplot(3,5,count); hist(ori_estimates(:,2), 0:180); hold on;
    plot([jj jj], get(gca, 'ylim'), 'r:');
    title([num2str(jj) '^{\circ}']); set(gca, 'xlim', [-1 181]);
    figure(f3); subplot(3,5,count); hist(ori_estimates(:,3), 0:180); hold on;
    plot([jj jj], get(gca, 'ylim'), 'r:');
    title([num2str(jj) '^{\circ}']); set(gca, 'xlim', [-1 181]);
    count = count+1;
end
