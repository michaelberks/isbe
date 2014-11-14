%%
% Do CLS detection on some nice mass regions


for ii = setdiff(17:53, [3 4 6 8 12 15 16 20]);% 24 34 36 40 51 52 53]
    
    %for each image
    i1 = double(imread(['C:\isbe\dev\background\images\mass_2\mass', zerostr(ii,3), '.bmp']));
    
    [gp p_sizes] = buildGpyr(i1, 5);
    gp = mb_change_pyramid_form(gp, p_sizes, 'g');
    
    CLS_results = cell(3,1);
    for jj = 1:3
        %Do CLS detection
        [CLS_results{jj}] = mb_cls_selection('ImageIn', gp{jj});
    end
    
    save(['C:\isbe\dev\background\misc\mass_cls_', zerostr(ii,3)], 'CLS_results');
    
%     cls_mr = CLS_results{1}.CLS | ...
%         imresize(CLS_results{2}.CLS, size(i1)) | ...
%             imresize(CLS_results{3}.CLS, size(i1));
%     
%     [r c] = find(cls_mr);    
%     figure('Name', ['mass', zerostr(ii,3)]); imagesc(i1);
%     colormap(gray(256)); axis image; hold on;
%     plot(c,r,'r.', 'MarkerSize', 4);
    
    clear i1 r c CLS_results cls_mr;
end
%%
% Do CLS detection on normal regions


for ii = 1:16 
    
    %for each image
    i1 = double(imread(['C:\isbe\dev\background\images\normal_2\normal', zerostr(ii,3), '.bmp']));
    
    [gp p_sizes] = buildGpyr(i1, 5);
    gp = mb_change_pyramid_form(gp, p_sizes, 'g');
    
    CLS_results = cell(3,1);
    for jj = 1:3
        %Do CLS detection
        [CLS_results{jj}] = mb_cls_selection('ImageIn', gp{jj});
    end
    
    save(['C:\isbe\dev\background\cls\normal_2\normal_cls_', zerostr(ii,3)], 'CLS_results');    
    clear i1 r c CLS_results cls_mr;
end
%%
for ii = 1:16
    load(['C:\isbe\dev\background\cls\normal_2\normal_cls_', zerostr(ii,3)]);
    i1 = double(imread(['C:\isbe\dev\background\images\normal_2\normal', zerostr(ii,3), '.bmp']));
    cls_mr = CLS_results{1}.CLS | ...
        imresize(CLS_results{2}.CLS, size(i1)) | ...
            imresize(CLS_results{3}.CLS, size(i1));

    figure; imagesc(i1); axis image; colormap(gray(256)); hold on;
    [r c] = find(cls_mr);
    plot(c, r, 'r.', 'MarkerSize', 4);
end

%%
for ii = [3 4 6 8 12 15 16 20]
    
    %for each image
    i1 = double(imread(['C:\isbe\dev\background\images\mass_2\mass', zerostr(ii,3), '.bmp']));
    
    load(['C:\isbe\dev\background\misc\mass_cls_', zerostr(ii,3)], 'CLS_results');
    
    [r1 c1] = find(CLS_results{1}.CLS);
    [r2 c2] = find(imresize(CLS_results{2}.CLS, size(i1)));
    [r3 c3] = find(imresize(CLS_results{3}.CLS, size(i1)));
       
    figure('Name', ['mass', zerostr(ii,3)]); imagesc(i1);
    colormap(gray(256)); axis image; hold on;
    plot(c3,r3,'g.', 'MarkerSize', 4);
    plot(c2,r2,'c.', 'MarkerSize', 4);
    plot(c1,r1,'r.', 'MarkerSize', 4);
    clear i1 r* c* CLS_results;
end
%%
figure;
subplot(1,2,1); imagesc(response1); axis image; colormap(gray(256));
subplot(1,2,2); imagesc(response2r); axis image; colormap(gray(256));
%%
figure;
subplot(1,2,1); imagesc(orientations1); axis image; colormap(gray(256));
subplot(1,2,2); imagesc(orientations2r); axis image; colormap(gray(256));
%%
figure;
subplot(1,2,1); imagesc(nms1); axis image; colormap(gray(256));
subplot(1,2,2); imagesc(nms2r); axis image; colormap(gray(256));
%%

display(['--Texture synthesis script started: ' datestr(now)]);
clear

for ii = 1:1
    pyr_args.TargetPyramid = u_load([mberksroot,...
        'background/g_pyramid/normal_2/normal', zerostr(ii, 3), '_pyramid.mat']);
    [rows cols] = size(pyr_args.TargetPyramid{1,1});
    row_centre = round(rows / 2);
    col_centre = round(cols / 2);
%     pyr_args.FilledImage = logical(ones(rows, cols)); %#ok
%     pyr_args.FilledImage(row_centre-64:row_centre+63, col_centre-64:col_centre+63) = 0;

    % Make a the biggest circular mask
    m = min([rows cols]);
    rad = floor((m - 128) / 2);
    [x y] = meshgrid(1:cols, 1:rows);
    pyr_args.FilledImage = (x - col_centre).^2 + (y - row_centre).^2 > rad.^2;
    clear x y m rad row_centre col_centre;

    pyr_args.ModelName = 'mass_2_g_20_model';

    pyr_args.ModelDir = [mberksroot, 'background/results/g_pyramid/'];
    pyr_args.CutOffLevel = 1;
    pyr_args.ConditionLevels = 1;
    pyr_args.SaveFile = [mberksroot, 'background/syn/g_pyramid/gp_circle', zerostr(ii, 3)];

    mb_gmm_g_pyr_synthesis(pyr_args);

    clear pyr_args;
end

clear
display(['--Texture synthesis script finished: ' datestr(now)]);
%%
for ii = 1:15
    load([mberksroot, 'background/syn/mass_2_2levels/conditioned_synthesis_circle_6_', zerostr(ii, 3)]);
    [p pp] = mb_change_pyramid_form(pyramid);
    synthesised_image = reconSFpyr(p, pp);
    save([mberksroot, 'background/syn/mass_2_2levels/conditioned_synthesis_circle_6_', zerostr(ii, 3)],...
        'synthesised_image', 'pyramid', 'cluster_image');
    figure; imagesc(synthesised_image); colormap(gray(256)); axis image;
end