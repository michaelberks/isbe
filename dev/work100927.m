sample_orientations1 = zeros(180,1);
sample_orientations2 = zeros(180,1);
for ii = 1:180
    cls = create_gauss_bar(3, 10, ii, 256, 256, 128, 128);
    [d d orientation_image local_ori] = monogenic_phase_cong(cls, 2, 8, 2, 0.65);
    sample_orientations1(ii) = orientation_image(126,126);
    sample_orientations2(ii) = 180*local_ori(126,126,2)/pi;
end
%%
for ii = 1:30
    
    load(['C:\isbe\asymmetry_project\data\synthetic_lines\lines512\image' zerostr(ii,3)]);
    orientation_image = u_load(['C:\isbe\asymmetry_project\data\synthetic_lines\lines512\results\monogenic_0p65_08_02_ori\probability_image' zerostr(ii,3)]);
    orientation_image(~(label==1)) = NaN;

    figure;
    subplot(1,2,1); imagesc(label_orientation); axis image; colormap([0 0 0; hsv(180)]); caxis([0 180]);
    subplot(1,2,2); imagesc(orientation_image); axis image; colormap([0 0 0; hsv(180)]); caxis([0 180]);
end
%%
for ii = 81:100
    
    prob_image1 = u_load(['C:\isbe\asymmetry_project\data\synthetic_lines\lines512\results\191877\probability_image' zerostr(ii,3)]);
    prob_image2 = u_load(['C:\isbe\asymmetry_project\data\synthetic_lines\lines512\results\192007\probability_image' zerostr(ii,3)]);

    figure;
    subplot(1,2,1); imagesc(prob_image1); axis image; caxis([0 1]);
    subplot(1,2,2); imagesc(prob_image2); axis image; caxis([0 1]);
end
%%
for ii = 2
    load(['C:\isbe\asymmetry_project\data\synthetic_lines\lines512\image' zerostr(ii,3)]);
    prob_image1 = u_load(['C:\isbe\asymmetry_project\data\synthetic_lines\lines512\results\191908\probability_image' zerostr(ii,3)]);
    prob_image2 = u_load(['C:\isbe\asymmetry_project\data\synthetic_lines\lines512\results\191906\probability_image' zerostr(ii,3)]);
    prob_image3 = u_load(['C:\isbe\asymmetry_project\data\synthetic_lines\lines512\results\192007\probability_image' zerostr(ii,3)]);
    
    [r1 a1(ii,1)] = calculate_roc_curve(prob_image1(:), label_centre(:),(-1:100)/100);
    [r2 a2(ii,1)] = calculate_roc_curve(prob_image2(:), label_centre(:),(-1:100)/100);
    [r3 a3(ii,1)] = calculate_roc_curve(prob_image3(:), label_centre(:),(-1:100)/100);

end
%%
figure; plot(r1(:,1), r1(:,2)); axis equal; axis([0 1 0 1]);
hold on;
plot(r2(:,1), r2(:,2), 'r'); axis equal; axis([0 1 0 1]);
plot(r3(:,1), r3(:,2), 'g'); axis equal; axis([0 1 0 1]);
plot(r1(:,1), r1(:,2), 'x'); axis equal; axis([0 1 0 1]);
plot(r2(:,1), r2(:,2), 'rx'); axis equal; axis([0 1 0 1]);
plot(r3(:,1), r3(:,2), 'gx'); axis equal; axis([0 1 0 1]);

tc = sum(label_centre(:));
fc = numel(label_centre)-tc;
for jj = 0:0.1:1
%     figure;
%     subplot(2,2,2); imagesc(prob_image1 > jj); axis image; caxis([0 1]);
%     subplot(2,2,3); imagesc(prob_image2 > jj); axis image; caxis([0 1]);
%     subplot(2,2,4); imagesc(prob_image3 > jj); axis image; caxis([0 1]);
%     subplot(2,2,1); imagesc(label_centre); axis image; caxis([0 1]);
    
    tp1 = sum((prob_image1(:) > jj) & label_centre(:)) / tc;
    fp1 = sum((prob_image1(:) > jj) & ~label_centre(:)) / fc;
    text(fp1, tp1, num2str(jj), 'color', 'b');
    
    tp2 = sum((prob_image2(:) > jj) & label_centre(:)) / tc;
    fp2 = sum((prob_image2(:) > jj) & ~label_centre(:)) / fc;
    text(fp2, tp2, num2str(jj), 'color', 'r');
    
    tp3 = sum((prob_image3(:) > jj) & label_centre(:)) / tc;
    fp3 = sum((prob_image3(:) > jj) & ~label_centre(:)) / fc;
    text(fp3, tp3, num2str(jj), 'color', 'g');
end
%%
for jj = 0:0.1:1
    figure;
    subplot(2,2,2); imagesc(prob_image1 > jj); axis image; caxis([0 1]);
    subplot(2,2,3); imagesc(prob_image2 > jj); axis image; caxis([0 1]);
    subplot(2,2,4); imagesc(prob_image3 > jj); axis image; caxis([0 1]);
    subplot(2,2,1); imagesc(label_centre); axis image; caxis([0 1]);
end
%%
s_args.num_levels = 1;
s_args.feature_shape = 'rect';
s_args.feature_type = 'all';
s_args.do_max = 0;
s_args.rotate = 0;
s_args.win_size = 3;
s_args.use_nag = 1;
classify_image(...
    'image_in', test_image,...
    'sampling_args', s_args,...
    'forest', [],...
    'decomp_type', 'dt',...
    'forest_type', 'isbe',...
    'use_probs', 0,...
    'mask', [],...
    'num_trees', [],...
    'max_size', 4);
%%
job = '192007';
for ii = 1:2
    s1 = load(['C:\isbe\asymmetry_project\data\synthetic_lines\lines512\image' zerostr(ii,3)]);
    s2 = load(['Z:\asymmetry_project\data\synthetic_lines\lines512\image' zerostr(ii,3)]);
    prob_image1 = u_load(['C:\isbe\asymmetry_project\data\synthetic_lines\lines512\results\' job '\probability_image' zerostr(ii,3)]);
    prob_image2 = u_load(['Z:\asymmetry_project\data\synthetic_lines\lines512\results\' job '\probability_image' zerostr(ii,3)]);
    [y1 x1] = find(s1.label_centre);
    [y2 x2] = find(s2.label_centre);
    display(sum(abs(prob_image1(:)-prob_image2(:))));
    
    figure; 
    subplot(1,2,1); imagesc(prob_image1 > 0.5); axis image; caxis([0 1]); hold on;
    plot(x1, y1, 'r.');
    plot(x2, y2, 'g.');
    subplot(1,2,2); imagesc(prob_image2 > 0.5); axis image; caxis([0 1]); hold on;
    plot(x1, y1, 'r.');
    plot(x2, y2, 'g.');

end
%%
mam_names = get_mammo_info(dir('Z:\asymmetry_project\data\mammograms\2004_screening\abnormals\mat\results\CC_01\*.mat'));
for ii = 1:length(mam_names)
    load(['Z:\asymmetry_project\data\mammograms\2004_screening\abnormals\mat\results\CC_01\' mam_names{ii} '_class.mat'])
    %load(['C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\mat\o04_' mam_names{ii} '.mat']);
    load(['C:\isbe\asymmetry_project\data\masks\2004_screening\abnormals\o04_' mam_names{ii} '_mask.mat']);
    load(['C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\mat\meta\' mam_names{ii} '_meta.mat']);
    
    if strcmpi(mam_names{ii}(4), 'R')
        probability_image = fliplr(probability_image);
    end
    
    figure; 
    subplot(1,2,1); imagesc(probability_image); axis image; caxis([0 1]);
    hold on;
    plot(meta_xy(:,1), meta_xy(:,2), 'r');
    subplot(1,2,2); hist(probability_image(mask), 101);
    display(['Mean classification probability of ' mam_names{ii} ' = ' num2str(mean(probability_image(mask)))]);
    clear probability_image mask meta_xy;
    saveas(gcf, ['M:\asymmetry_project\Weekly presentations\figures\abnormal' zerostr(ii,2) '.tif']);
    close(gcf);    
    
end
%
mam_names = get_mammo_info(dir('Z:\asymmetry_project\data\mammograms\2004_screening\normals\mat\results\CC_01\*.mat'));
for ii = 1:length(mam_names)
    load(['Z:\asymmetry_project\data\mammograms\2004_screening\normals\mat\results\CC_01\' mam_names{ii} '_class.mat'])
    %load(['C:\isbe\asymmetry_project\data\mammograms\2004_screening\normals\mat\' mam_names{ii} '.mat']);
    load(['C:\isbe\asymmetry_project\data\masks\2004_screening\normals\' mam_names{ii} '_mask.mat']);
    
    if strcmpi(mam_names{ii}(4), 'R')
        probability_image = fliplr(probability_image);
    end
    
    figure; 
    subplot(1,2,1); imagesc(probability_image); axis image; caxis([0 1]);
    subplot(1,2,2); hist(probability_image(mask), 101);
    display(['Mean classification probability of ' mam_names{ii} ' = ' num2str(mean(probability_image(mask)))]);
    clear probability_image mask meta_xy;
    saveas(gcf, ['M:\asymmetry_project\Weekly presentations\figures\normal' zerostr(ii,2) '.tif']);
    close(gcf);
end