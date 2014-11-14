%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Testing how image-pair RF classification performs given rotated copies of
% the same region
%--------------------------------------------------------------------------
%%
% First test: do two-image RF classification on the same image
norm_list = dir('C:\isbe\dev\background\images\normal512\*.bmp');

i1 = double(imread(['C:\isbe\dev\background\images\normal512\' norm_list(1).name]));

args.image1 = i1;
args.image2 = i1;
args.num_train = 10000;
args.num_levels = 5;
args.d = 5;
args.tree_dir = 'C:\isbe\dev\misc';
args.do_test1 = 1;
args.do_test2 = 0;
args.use_probs = 0;
args.do_max = 0;
args.win_size = 1;
args.max_test_size = 64;
args.n_trees = 100;
args.save_path = 'C:\isbe\asymmetry_project\misc\same_image_rf.mat';
clear i1;
%
[random_forest] =...
    mb_random_forest_two_images(args);
%
im1_votes1 = random_forest.image1_votes1 ./ random_forest.image1_total_votes;
im1_votes1(~random_forest.image1_total_votes) = 0;

figure; imagesc(im1_votes1); axis image; colorbar;

%%
% Now do two image RF classification on rotated copies of the same region -
% note that the script below has been submitted to hydra to run
im = double(imread([asymmetryroot 'data/contralateral_rfs/misc/o04_010RCC_1024_3797_3365.bmp']));
xx = repmat(1:512, 512, 1);
mask = (xx-256).^2 + (xx' - 256).^2 < 256^2;
for theta = 9*(0:19)

    i1 = imrotate(im, theta, 'crop');
    i1 = i1(257:768, 257:768);

    i2 = imrotate(im, -theta, 'crop');
    i2 = i2(257:768, 257:768);

    args.image1 = i1;
    args.image2 = i2;
    args.mask1 = mask;
    args.mask2 = mask;
    args.num_train = 10000;
    args.num_levels = 5;
    args.d = 20;
    args.tree_dir = [asymmetryroot 'data/misc/'];
    args.do_test1 = 1;
    args.do_test2 = 1;
    args.use_probs = 0;
    args.do_max = 0;
    args.win_size = 3;
    args.max_test_size = 64;
    args.n_trees = 100;
    args.save_path = [asymmetryroot 'data/misc/same_image_rf' zerostr(2*theta, 3) '.mat'];
    clear i1 i2;
    %
    mb_random_forest_two_images(args);

end
%%
im = double(imread([asymmetryroot 'data/contralateral_rfs/misc/o04_010RCC_1024_3797_3365.bmp']));
code = 'W1L5M0';
mkdir('C:\isbe\asymmetry_project\data\contralateral_rfs\misc\figures\ilp');
for theta = [3 4 7]%9*(1:5)
    
    i1 = imrotate(im, theta, 'crop');
    i1 = i1(257:768, 257:768);

    i2 = imrotate(im, -theta, 'crop');
    i2 = i2(257:768, 257:768);

    load([asymmetryroot 'data/contralateral_rfs/misc/ilp/same_image_rf_' code 't' zerostr(2*theta, 3) '.mat']);

    im1_votes1 = random_forest.image1_votes1 ./ random_forest.image1_total_votes;
    im1_votes1(~random_forest.image1_total_votes) = NaN;

    im2_votes1 = random_forest.image2_votes1 ./ random_forest.image2_total_votes;
    im2_votes1(~random_forest.image2_total_votes) = NaN;
    
    hist_min = min([min(im1_votes1(:)) min(im2_votes1(:))]);
    hist_max = max([max(im1_votes1(:)) max(im2_votes1(:))]);
    hist1 = hist(im1_votes1(:), hist_min:0.01:hist_max);
    hist2 = hist(im2_votes1(:), hist_min:0.01:hist_max);
    
    [roc_pts auc] = calculate_roc_curve(...
        [im1_votes1(random_forest.image1_total_votes>0); im2_votes1(random_forest.image2_total_votes>0)],...
        [true(sum(random_forest.image1_total_votes(:)>0),1); false(sum(random_forest.image2_total_votes(:)>0),1)],...
        (-1:101)/100);
    
    display( [...
        'Angle diff = ', num2str(2*theta)...
        ' mean = ' num2str(naNmean(im1_votes1(:)),3) ', ' num2str(naNmean(im2_votes1(:)),3) ', ' num2str(naNmean(im1_votes1(:))-naNmean(im2_votes1(:)),3) '; '...
        ' auc = ' num2str(auc,3)] );
    
    if 0%ismember(theta, 9*[1 5])
        code = 'ilp';
        write_im_from_colormap(i1, ['C:\isbe\asymmetry_project\data\contralateral_rfs\misc\figures\' code '\im1_' zerostr(2*theta, 3) '.bmp'], gray(256));
        write_im_from_colormap(i2, ['C:\isbe\asymmetry_project\data\contralateral_rfs\misc\figures\' code '\im2_' zerostr(2*theta, 3) '.bmp'], gray(256));
        write_im_from_colormap(im1_votes1, ['C:\isbe\asymmetry_project\data\contralateral_rfs\misc\figures\' code '\im1v1_' zerostr(2*theta, 3) '.bmp'], jet(256));
        write_im_from_colormap(im2_votes1, ['C:\isbe\asymmetry_project\data\contralateral_rfs\misc\figures\' code '\im2v1_' zerostr(2*theta, 3) '.bmp'], jet(256));

        f1 = figure('windowstyle', 'normal', 'units', 'pixels', 'position', [50 50 800 800]);
        bar(hist_min:0.01:hist_max, [hist1' hist2']); set(gca, 'XLim', [hist_min hist_max]); colormap([0 0 1; 1 0 0]);
        title(['RF classification - histogram of probabilities \theta = \pm' num2str(theta)]);
        saveas(f1, ['C:\isbe\asymmetry_project\data\contralateral_rfs\misc\figures\' code '\hist_' zerostr(2*theta, 3) '.bmp']);
        close(f1);
        
        f1 = figure('windowstyle', 'normal', 'units', 'pixels', 'position', [50 50 800 800]); 
        plot(roc_pts(:,1), roc_pts(:,2), 'linewidth', 2); axis equal; axis([0 1 0 1]);
        legend(['AUC = ' num2str(auc,3)], 'location', 'southeast');
        title(['ROC curve for rotation \theta = \pm' num2str(theta)]);
        xlabel('FPF');
        ylabel('TPF');
        saveas(f1, ['C:\isbe\asymmetry_project\data\contralateral_rfs\misc\figures\' code '\roc_' zerostr(2*theta, 3) '.bmp']);
        close(f1);
        code = 'W1L5M0';
    end
    
%     figure; 
%     subplot(1,2,1); imagesc(i1); axis image;
%     subplot(1,2,2); imagesc(i2); axis image; colormap(gray(256));
    figure; 
    subplot(1,2,1); imagesc(im1_votes1); axis image;
    subplot(1,2,2); imagesc(im2_votes1); axis image; colormap(jet(256));
    figure; 
    subplot(1,2,1); bar(hist_min:0.01:hist_max, hist1, 'b'); set(gca, 'XLim', [hist_min hist_max]);
    subplot(1,2,2); bar(hist_min:0.01:hist_max, hist2, 'b'); set(gca, 'XLim', [hist_min hist_max]);
end
%%
mkdir('C:\isbe\asymmetry_project\data\contralateral_rfs\misc\figures\ilp');
for offset = [1 2 3]%9*(1:5)

    load([asymmetryroot 'data/contralateral_rfs/misc/offset/off_image_rf_W1L5M0t' zerostr(offset, 2) '.mat']);

    im1_votes1 = random_forest.image1_votes1 ./ random_forest.image1_total_votes;
    im1_votes1(~random_forest.image1_total_votes) = NaN;

    im2_votes1 = random_forest.image2_votes1 ./ random_forest.image2_total_votes;
    im2_votes1(~random_forest.image2_total_votes) = NaN;
    
    hist_min = min([min(im1_votes1(:)) min(im2_votes1(:))]);
    hist_max = max([max(im1_votes1(:)) max(im2_votes1(:))]);
    hist1 = hist(im1_votes1(:), hist_min:0.01:hist_max);
    hist2 = hist(im2_votes1(:), hist_min:0.01:hist_max);
    
    [roc_pts auc] = calculate_roc_curve(...
        [im1_votes1(random_forest.image1_total_votes>0); im2_votes1(random_forest.image2_total_votes>0)],...
        [true(sum(random_forest.image1_total_votes(:)>0),1); false(sum(random_forest.image2_total_votes(:)>0),1)],...
        (-1:101)/100);
    
    display( [...
        'Offset = ', num2str(offset)...
        ' mean = ' num2str(naNmean(im1_votes1(:)),3) ', ' num2str(naNmean(im2_votes1(:)),3) ', ' num2str(naNmean(im1_votes1(:))-naNmean(im2_votes1(:)),3) '; '...
        ' auc = ' num2str(auc,3)] );
    
    if 0%ismember(theta, 9*[1 5])
        code = 'ilp';
        write_im_from_colormap(i1, ['C:\isbe\asymmetry_project\data\contralateral_rfs\misc\figures\' code '\im1_' zerostr(2*theta, 3) '.bmp'], gray(256));
        write_im_from_colormap(i2, ['C:\isbe\asymmetry_project\data\contralateral_rfs\misc\figures\' code '\im2_' zerostr(2*theta, 3) '.bmp'], gray(256));
        write_im_from_colormap(im1_votes1, ['C:\isbe\asymmetry_project\data\contralateral_rfs\misc\figures\' code '\im1v1_' zerostr(2*theta, 3) '.bmp'], jet(256));
        write_im_from_colormap(im2_votes1, ['C:\isbe\asymmetry_project\data\contralateral_rfs\misc\figures\' code '\im2v1_' zerostr(2*theta, 3) '.bmp'], jet(256));

        f1 = figure('windowstyle', 'normal', 'units', 'pixels', 'position', [50 50 800 800]);
        bar(hist_min:0.01:hist_max, [hist1' hist2']); set(gca, 'XLim', [hist_min hist_max]); colormap([0 0 1; 1 0 0]);
        title(['RF classification - histogram of probabilities \theta = \pm' num2str(theta)]);
        saveas(f1, ['C:\isbe\asymmetry_project\data\contralateral_rfs\misc\figures\' code '\hist_' zerostr(2*theta, 3) '.bmp']);
        close(f1);
        
        f1 = figure('windowstyle', 'normal', 'units', 'pixels', 'position', [50 50 800 800]); 
        plot(roc_pts(:,1), roc_pts(:,2), 'linewidth', 2); axis equal; axis([0 1 0 1]);
        legend(['AUC = ' num2str(auc,3)], 'location', 'southeast');
        title(['ROC curve for rotation \theta = \pm' num2str(theta)]);
        xlabel('FPF');
        ylabel('TPF');
        saveas(f1, ['C:\isbe\asymmetry_project\data\contralateral_rfs\misc\figures\' code '\roc_' zerostr(2*theta, 3) '.bmp']);
        close(f1);
        code = 'W1L5M0';
    end
    
%     figure; 
%     subplot(1,2,1); imagesc(i1); axis image;
%     subplot(1,2,2); imagesc(i2); axis image; colormap(gray(256));
    figure; 
    subplot(1,2,1); imagesc(im1_votes1); axis image;
    subplot(1,2,2); imagesc(im2_votes1); axis image; colormap(jet(256));
    figure; 
    subplot(1,2,1); bar(hist_min:0.01:hist_max, hist1, 'b'); set(gca, 'XLim', [hist_min hist_max]);
    subplot(1,2,2); bar(hist_min:0.01:hist_max, hist2, 'b'); set(gca, 'XLim', [hist_min hist_max]);
end