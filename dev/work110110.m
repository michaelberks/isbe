test_dir = 'C:\isbe\asymmetry_project\data\synthetic_lines\real512\';
prob_dir = 'C:\isbe\asymmetry_project\data\synthetic_lines\real512\results\';

%%
%Compute line and orientation maps for analytic g2d methods
param_dir = 'g2d_scales_16';

mkdir([prob_dir param_dir '\lines\']);
mkdir([prob_dir param_dir '\orientations\']);

for ii = 1:100
    display(['Processing image ', num2str(ii) ' of 100']);
    %Load in test image
    test_im = u_load([test_dir '\image' zerostr(ii,3) '.mat']);
    
    [line_orientation, line_map] = karssemeijer_line_detection(...
        test_im,...
        'line_scales', [1 2 4 8],...
        'grad_scale', 10,...
        'grad_ori_thresh', pi/6,...
        'grad_strength_thresh', 25,...
        'line_strength_thresh', 0,...
        'binary_map', 0);
    
    line_orientation = 180*line_orientation/pi;
    if ii < 11
        figure; colormap(gray(256));
        subplot(1,2,1); imagesc(test_im); axis image;
        subplot(1,2,2); imagesc(line_map); axis image;
    end
        
    save([prob_dir param_dir '\lines\image' zerostr(ii,3) '_class.mat'], 'line_map');
    save([prob_dir param_dir '\orientations\image' zerostr(ii,3) '_class.mat'], 'line_orientation');

end
%%
param_dir = 'g2d_scales_16\orientations\';
[orientation_errors] = compute_image_orientation_errors(...
    test_dir, [prob_dir param_dir], 'centre_line');
display(mean(abs(orientation_errors)));

param_dir = '233902';
[orientation_errors] = compute_image_orientation_errors(...
    test_dir, [prob_dir param_dir], 'centre_line');
display(mean(abs(orientation_errors)));

param_dir = '233908';
[orientation_errors] = compute_image_orientation_errors(...
    test_dir, [prob_dir param_dir], 'centre_line');
display(mean(abs(orientation_errors)));
%%
param_dirs = {'233902', '233908', 'g2d_scales_16\orientations\'};
titles = {'RF/DT-CWT, real backgrounds', 'RF/DT-CWT, synthetic backgrounds', 'Gaussian 2nd derivatives'};

spacing = 2;
num_angles = 36;
ang_res = 180 / num_angles;
arrow_colors = hsv(num_angles);

for ii = 1:10
    test_im = u_load([test_dir 'image' zerostr(ii,3) '.mat']);
    label = load([test_dir 'labels\label' zerostr(ii,3) '.mat']);
    
    a = zeros(4,1);
    figure;
    
    for kk = 1:4
        if kk == 4
            orientation_map = label.label_orientation;
            a(kk) = subplot(2,2,kk); imagesc(label.label); axis image; colormap(gray(256)); hold on;
        else
            orientation_map = load_uint8([prob_dir param_dirs{kk} '\image' zerostr(ii,3) '_class.mat']);
            a(kk) = subplot(2,2,kk); imagesc(test_im); axis image; colormap(gray(256)); hold on;
            
            error_mask = label.label == 1; 
            ori_errors = mb_mod(orientation_map(error_mask) - label.label_orientation(error_mask), 180);
            title(titles{kk});
            xlabel(['Mean orientation error = ' num2str(mean(abs(ori_errors)))]);
        end
    
        
        
        mask = label.label > 0;
        [r c] = size(mask);

        spacing_mask = false(r, c);
        spacing_mask(1:spacing:r, 1:spacing:c) = true;
        mask = mask & spacing_mask;

        for jj = 0:num_angles
            theta = (jj - 0.5)*ang_res;

            %Get mask of pixels that have orientation within theta range
            angle_mask = mask & ...
                 (orientation_map > theta - 0.5*ang_res) &...
                 (orientation_map <= theta + 0.5*ang_res);

            [y x] = find(angle_mask);
            u = cos(pi*orientation_map(angle_mask)/180);
            v = -sin(pi*orientation_map(angle_mask)/180);

            quiver(a(kk), x, y, 4*u, 4*v, 0, 'color', arrow_colors(mod(jj, num_angles)+1,:));
        end
    end
    linkaxes(a);
    zoom on;
end