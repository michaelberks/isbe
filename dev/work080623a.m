% Work 23rd June 2008 - First m-file on my new computer!!

[gp1 p_sizes] = buildGpyr(i1, 5);
    gp1 = mb_change_pyramid_form(gp1, p_sizes, 'g');
    
for level = 1:5
    figure; imagesc(gp1{level}); axis image; colormap(gray(256));
end

[CLS_result.GaborResponse CLS_result.GaborOrientation] = ...
    gabor_filter('ImageIn', gp1{4}, 'Tau', 2, 'Len', 5);

%%
GaborFilterArgs.Tau = 2;
GaborFilterArgs.Len = 5;
%%
[CLS_result] = mb_cls_selection('ImageIn', gp1{4}, 'GaborFilterArgs', GaborFilterArgs, 'MinCLSLength', 8);
figure; imagesc(CLS_result.CLS); axis image; colormap(gray(256));
%%
im_list = dir('C:\isbe\dev\background\images\normal512\o04_*');
dt_list = dir('C:\isbe\dev\background\dual_tree\normal512\*tree*');
%%

GaborFilterArgs.NumAngles = 12;
GaborFilterArgs.Tau = 4;
GaborFilterArgs.Len = 16;
%%
orientations = [Inf -Inf; 0 pi; pi 2*pi; 2*pi 3*pi; -3*pi -2*pi; -2*pi -pi; -pi 0] / 6;
my_map = [0 0 0; jet(7)];
%%
%[lengths, tau] = meshgrid([4, 8], [2, 4]);
for ii = 1:10
    
    %for each image
    i1 = double(imread(['C:\isbe\dev\background\images\normal512\', im_list(ii).name]));
    
    [gp p_sizes] = buildGpyr(i1, 5);
    gp = mb_change_pyramid_form(gp, p_sizes, 'g');    
    
    dt = u_load(['C:\isbe\dev\background\dual_tree\normal512\', dt_list(jj).name]);
    
    CLS_result = cell(1);
    
    for level = 2:2
        
        [CLS_result{1}] = mb_cls_selection('ImageIn', gp{level+1}, 'GaborFilterArgs', GaborFilterArgs, 'MinCLSLength', 10*(5-level), 'Connectivity', 0);
        
        [CLS_map] =  mb_cls_map_orientations(CLS_result, orientations);
        figure;
        for ori = 1:6
            [y_pts x_pts] = find(CLS_map{1} == ori+1);
            subplot(2,3,ori); imagesc(real(dt{level}(:,:,ori))); axis image; colormap(gray(256)); hold on;
            plot(x_pts, y_pts, 'r.', 'MarkerSize', 4);
            %title(['Tau = ', num2str(tau(jj)), ' Length = ', num2str(lengths(jj))]);
            title(['Level = ', num2str(level), ' Ori = ', num2str(ori)]);
        end
    end
        
    clear i1 gp CLS_result;
end
%%
im_list = dir('C:\isbe\mammograms\new_CAD\BMP_2004\*.bmp');
%
mkdir('C:\isbe\mammograms\new_CAD\BMP_2004_half\');
for ii = 1:length(im_list);
    
    display(['resizing mammogram ', num2str(ii)]);
    i1 = imread(['C:\isbe\mammograms\new_CAD\BMP_2004\', im_list(ii).name]);
    mammogram = imresize(i1, 0.5, 'bilinear');
    clear i1;
    save(['C:\isbe\mammograms\new_CAD\BMP_2004_half\', im_list(ii).name(1:end-4)], 'mammogram');
    clear mammogram;
end
%%
m_list = dir('C:\isbe\mammograms\new_CAD\BMP_2004_half\*.mat');
mask_list = dir('C:\isbe\mammograms\new_CAD\BMP_2004_masks\*.mat');
%%
test_im_width = ones(512)*120;
for ii = 1:7
    test_im_width(207:306, ii*64:(ii*64)+ii-1) = 200;
end
%%
test_im_length = ones(512)*120;
for ii = 1:7
    test_im_length(256-(ii*5)+1:256+(ii*5), ii*64:(ii*64)+1) = 200;
end    


%%
sum(nms4l16(:) & (test_im_length(:) == 200)) /...
    ( sum(nms4l16(:) & (test_im_length(:) == 200)) + ...
     ~sum(nms4l16(:) & (test_im_length(:) == 200)) + ...
      sum(nms4l16(:) & ~(test_im_length(:) == 200)) );