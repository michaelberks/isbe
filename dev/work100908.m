generate_synthetic_curve_images('num_images', 100, 'save_dir', 'curves512/', 'plot', 0);
generate_synthetic_curve_images('num_images', 100, 'save_dir', 'lines512/', 'line_type', 'sin', 'plot', 0);
%%
bg = u_load('C:\isbe\asymmetry_project\data\synthetic_backgrounds\smooth512\train\bg00500.mat');
cls = create_gauss_bar(3, 10, 0, 512, 512, 256.5, 256.5);
cls = cls+bg;

mps = zeros(360, 12*7*4);
mps_aligned = zeros(360, 12*7*4);

for ii = 1:360
    
        %angled_cls = imrotate(cls, ii,'crop', 'bicubic');
        %angled_cls = angled_cls(129:384, 129:384);
        angled_cls = create_gauss_bar(2, 10, ii, 256, 256, 128.5, 128.5);
        dt = dtwavexfm2b(angled_cls, 4);
        
        mps(ii,:) = sample_dt_data(dt, 128.5, 128.5, 'feature_shape', 'clock', 'feature_type', 'conj', 'rotate', 0);
        mps_aligned(ii,:) = sample_dt_data(dt, 128.5, 128.5, 'feature_shape', 'clock', 'feature_type', 'conj', 'rotate', 1);

end
%
%%
jump = 30;

max_mag = max(abs(mps(:)));
letters = 'MABCDEF';
for lev = 3
    for pt = 1:7
        figure('name', ['DT magnitudes at centre of rotation for point ' letters(pt)]);
        for ori = 1:6
            subplot(2,3,ori); hold all;
            
            col = (lev-1)*42 + (ori-1)*7 + pt;
            plot(1:jump:360, mps(1:jump:360,col)); axis([1 360 0 max_mag]);
            plot(1:jump:360, mps_aligned(1:jump:360,col), 'r--'); axis([1 360 0 max_mag]);

            xlabel('Angle (degrees)');
            ylabel('DT-CWT magnitude');
            title(['Band ' num2str(ori)]);

%             if ii == 1 || ii == 6
%                 loc = 'southeast';
%             else
%                 loc = 'northeast';
%             end
%             legend({...
%                 ['Interpolated, std = ' num2str(std(mags_interp(:,ii)))],...
%                 ['Magnitude aligned, std = ' num2str(std(mags_aligned(:,ii)))],...
%                 ['Complex aligned, std = ' num2str(std(abs(comps_aligned(:,ii))))],...
%                 ['Phase and mag aligned, std = ', num2str(std(mps_aligned(:,ii)))]},...
%                 'location' , loc);
        end

        figure('name', ['DT phases at centre of rotation for point ' letters(pt)]);
        for ori = 1:6
            subplot(2,3,ori); hold all;
            
            col = (lev+3)*42 + (ori-1)*7 + pt;
            
            plot(1:jump:360, mps(1:jump:360,col)); axis([1 360 -pi pi]);
            plot(1:jump:360, mps_aligned(1:jump:360,col), 'r--'); axis([1 360 -pi pi]);

            xlabel('Angle (degrees)');
            ylabel('DT-CWT phase');
            title(['Band ' num2str(ori)]);

%             legend({...
%                 ['Interpolated, std = ' num2str(std(phases_interp(:,ii)))],...
%                 ['Phase aligned, std = ' num2str(std(phases_aligned(:,ii)))],...
%                 ['Complex aligned, std = ' num2str(std(angle(comps_aligned(:,ii))))],...
%                 ['Phase and mag aligned, std = ', num2str(std(mps_aligned(:,ii+6)))]},...
%                 'location' , 'southeast');
        end
    end
end
%%
bg = u_load('C:\isbe\asymmetry_project\data\synthetic_backgrounds\smooth512\train\bg00500.mat');
cls = create_gauss_bar(3, 10, 0, 512, 512, 256.5, 256.5);
cls = cls+bg;

mps = zeros(360, 12*4);
mps_aligned = zeros(360, 12*4);

for ii = 1:360
    
        %angled_cls = imrotate(cls, ii,'crop', 'bicubic');
        %angled_cls = angled_cls(129:384, 129:384);
        angled_cls = create_gauss_bar(2, 10, ii, 256, 256, 128.5, 128.5);
        dt = dtwavexfm2b(angled_cls, 4);
        
        mps(ii,:) = sample_dt_data(dt, 128.5, 128.5, 'feature_shape', 'rect', 'feature_type', 'conj', 'rotate', 0, 'win_size', 1);
        mps_aligned(ii,:) = sample_dt_data(dt, 128.5, 128.5, 'feature_shape', 'rect', 'feature_type', 'conj', 'rotate', 1, 'win_size', 1);

end
%
%
jump = 1;

max_mag = max(abs(mps(:)));
for lev = 1:4
    figure('name', ['DT magnitudes at centre of rotation for level ' num2str(lev)]);
    for ori = 1:6
        subplot(2,3,ori); hold all;

        col = (lev-1)*6 + ori;
        plot(1:jump:360, mps(1:jump:360,col)); axis([1 360 0 max_mag]);
        plot(1:jump:360, mps_aligned(1:jump:360,col)); axis([1 360 0 max_mag]);

        xlabel('Angle (degrees)');
        ylabel('DT-CWT magnitude');
        title(['Band ' num2str(ori)]);

%             if ii == 1 || ii == 6
%                 loc = 'southeast';
%             else
%                 loc = 'northeast';
%             end
%             legend({...
%                 ['Interpolated, std = ' num2str(std(mags_interp(:,ii)))],...
%                 ['Magnitude aligned, std = ' num2str(std(mags_aligned(:,ii)))],...
%                 ['Complex aligned, std = ' num2str(std(abs(comps_aligned(:,ii))))],...
%                 ['Phase and mag aligned, std = ', num2str(std(mps_aligned(:,ii)))]},...
%                 'location' , loc);
    end

    figure('name', ['DT phases at centre of rotation for level ' num2str(lev)]);
    for ori = 1:6
        subplot(2,3,ori); hold all;

        col = (lev+3)*6 + ori;

        plot(1:jump:360, mps(1:jump:360,col)); axis([1 360 -pi pi]);
        plot(1:jump:360, mps_aligned(1:jump:360,col)); axis([1 360 -pi pi]);

        xlabel('Angle (degrees)');
        ylabel('DT-CWT phase');
        title(['Band ' num2str(ori)]);

%             legend({...
%                 ['Interpolated, std = ' num2str(std(phases_interp(:,ii)))],...
%                 ['Phase aligned, std = ' num2str(std(phases_aligned(:,ii)))],...
%                 ['Complex aligned, std = ' num2str(std(angle(comps_aligned(:,ii))))],...
%                 ['Phase and mag aligned, std = ', num2str(std(mps_aligned(:,ii+6)))]},...
%                 'location' , 'southeast');
    end
end
%%
for ii = 0:5;
    figure; hold all;
    plot(1:360, mps_aligned(15,1+ii) ./ mps_aligned(:,1+ii), 'linewidth', 2); axis([0 360 0 2]);
    plot(1:360, mps_aligned(15,7+ii) ./ mps_aligned(:,7+ii), 'linewidth', 2); axis([0 360 0 2]);
    plot(1:360, mps_aligned(15,13+ii) ./ mps_aligned(:,13+ii), '--', 'linewidth', 2); axis([0 360 0 2]);
    plot(1:360, mps_aligned(15,19+ii) ./ mps_aligned(:,19+ii), '--', 'linewidth', 2); axis([0 360 0 2]);
end
%%
o_list = dir('Z:\asymmetry_project\data\synthetic_lines\curves512\results\191878\*.mat');
for ii = 1:30%length(o_list)
    load(['C:\isbe\asymmetry_project\data\synthetic_lines\curves512\', o_list(ii).name(end-11:end)]);
    reg_orientation = u_load(['Z:\asymmetry_project\data\synthetic_lines\curves512\results\191878\' o_list(ii).name]);
    reg_orientation(~(label==1)) = NaN;
    
    figure;
    subplot(1,2,1); imagesc(label_orientation); axis image; colormap([0 0 0; hsv(255)]); caxis([0 180]);
    subplot(1,2,2); imagesc(reg_orientation); axis image; colormap([0 0 0; hsv(255)]); caxis([0 180]);
    
end
%%
generate_pixel_training_data(...
    'num_samples', 2e3,...
    'bg_dir', [asymmetryroot, 'data/synthetic_backgrounds/smooth512/train/'], ...
    'num_bgs', [],...
    'bg_stem', [],...
    'bg_zeros', [],...
    'detection_type', 'detection',...
    'win_size', 15,...
    'bg_ratio', 1,...
    'pts_per_image', 400,...
    'width_range', [4 16],...
    'orientation_range', [0 360],...
    'contrast_range', [4 16],...
    'decay_rate', 4,...
    'line_type', 'sin',...
    'save_path', [], ...
    'plot', 0);
%%
[roc_p0, auc_p0] = compute_roc_image_set_lines('C:\isbe\asymmetry_project\data\synthetic_lines\lines512\', 'C:\isbe\asymmetry_project\data\synthetic_lines\lines512\results\191874\', 'centre_line');
[roc_p1, auc_p1] = compute_roc_image_set_lines('C:\isbe\asymmetry_project\data\synthetic_lines\lines512\', 'C:\isbe\asymmetry_project\data\synthetic_lines\lines512\results\191875\', 'centre_line');
[roc_m0, auc_m0] = compute_roc_image_set_lines('C:\isbe\asymmetry_project\data\synthetic_lines\lines512\', 'C:\isbe\asymmetry_project\data\synthetic_lines\lines512\results\191876\', 'centre_line');
[roc_m1, auc_m1] = compute_roc_image_set_lines('C:\isbe\asymmetry_project\data\synthetic_lines\lines512\', 'C:\isbe\asymmetry_project\data\synthetic_lines\lines512\results\191877\', 'centre_line');
