%%
% What are we trying to do here? Perform CLS detection on whole mammograms.
% However in lower levels artefacts are caused in the valid breast region
% due to filtering the non-valid region (e.g. the near white edges of the
% mammogram). Let's black the se regions out first;

%load mammogram and mask of breast region
m_list = dir('C:\isbe\mammograms\new_CAD\BMP_2004_half\*.mat');
mammogram = double(u_load(['C:\isbe\mammograms\new_CAD\BMP_2004_half\', m_list(3).name]));
mask = u_load(['C:\isbe\mammograms\new_CAD\BMP_2004_masks\', m_list(3).name(1:end-4), '_mask']);

%build Gaussian Pyramid of mammogram
[gp p_sizes] = buildGpyr(mammogram, 5);
gp = mb_change_pyramid_form(gp, p_sizes, 'g');

%Resize mask for each level of the pyramid and set non-breast region to
%zero
masks = cell(5,1);
for level = 1:5
    masks{level} = imresize(~mask, size(gp{level}));
    %gp{level}(masks{level}) = 0;
end
clear mask;

%Perform CLS detection on each level of the pyramid and plot result
cls = cell(5,1);
GaborFilterArgs.Normalise = 1;
%%
GaborFilterArgs.HighPassFilter = true;

%-imresize(gp{level+1}, size(gp{level}))

for level = 2:5
    cls{level} = mb_cls_selection('ImageIn', gp{level},...
        'GaborFilterArgs', GaborFilterArgs,...
        'MinLength', 12, 'Connectivity', 0, 'Thin', 1,...
        'IgnoreMap', [], 'GaborThreshold', 0,...
        'GradientThreshold', 0, 'GradientAlignment', pi/12,...
        'NMS', true, 'Debug', false);
    
    [y_pts x_pts] = find(cls{level}.CLS > 0);
    
    figure; 
    %subplot(1,2,1); imagesc(gp{level}); axis image; colormap(gray(256));
    imagesc(gp{level}); axis image; colormap(gray(256));
    hold on; plot(x_pts, y_pts, 'r.', 'MarkerSize', 4);
    
%     temp = cls{level}.GaborResponse;
%     temp(masks{level}) = 0;
%     subplot(1,2,2); imagesc(temp); axis image; colormap(gray(256));
%     hold on; plot(x_pts, y_pts, 'r.', 'MarkerSize', 4);
end
%%
% But now we get a really strong artefact due the sharp edge of the breast
% region! Maybe we need to try the gradient discard method used by
% Rangayyan
%Build Gaussian Derivative filter if not supplied

%Make Gaussian filter
sigma = 5;
width = 21; %Should probably give the user the option of choosing this

ssq = sigma^2;
[x,y] = meshgrid(-width:width,-width:width);
GaussDeriv = -x.*exp(-(x.*x+y.*y)/(2*ssq))/(pi*ssq);

for level = 3:5

    %Calulate image gradient and normal direction to gradient
    grad_x = imfilter(gp{level},GaussDeriv, 'conv','replicate');
    grad_y = imfilter(gp{level}, GaussDeriv', 'conv','replicate');

    normal_to_gradient = atan(-grad_x ./ grad_y);
    gradient_mag = grad_x.^2 + grad_y.^2;
    
    a = 2;%(6 - level)*2;
    
%     figure; 
%     
%     subplot(1,2,1); imagesc(gp{level}); axis image; colormap(gray(256));
%     hold on;
%     %quiver(1:a:size(gp{level},2), 1:a:size(gp{level},1), grad_x(1:a:end, 1:a:end), grad_y(1:a:end, 1:a:end));
%     quiver(1:a:size(gp{level},2), 1:a:size(gp{level},1), cos(normal_to_gradient(1:a:end, 1:a:end)), sin(normal_to_gradient(1:a:end, 1:a:end)));
%     a1 = gca;
%     
%     subplot(1,2,2); imagesc(normal_to_gradient + cls{level}.GaborOrientation); axis image; colormap(hsv(256));
%     quiver(1:a:size(gp{level},2), 1:a:size(gp{level},1), cos(-cls{level}.GaborOrientation(1:a:end, 1:a:end)), sin(-cls{level}.GaborOrientation(1:a:end, 1:a:end)), 'r');
%     hold on;
%     
%     a2 = gca;
%     linkaxes([a1, a2]);
    
    figure; imagesc(gp{level}); axis image; colormap(gray(256));
    hold on;
    quiver(1:a:size(gp{level},2), 1:a:size(gp{level},1), cos(normal_to_gradient(1:a:end, 1:a:end)), sin(normal_to_gradient(1:a:end, 1:a:end)));
    quiver(1:a:size(gp{level},2), 1:a:size(gp{level},1), cos(-cls{level}.GaborOrientation(1:a:end, 1:a:end)), sin(-cls{level}.GaborOrientation(1:a:end, 1:a:end)), 'r');
    
    figure; imagesc(normal_to_gradient + cls{level}.GaborOrientation); axis image; colormap(hsv(256));
    hold on;
    quiver(1:a:size(gp{level},2), 1:a:size(gp{level},1), cos(normal_to_gradient(1:a:end, 1:a:end)), sin(normal_to_gradient(1:a:end, 1:a:end)));
    quiver(1:a:size(gp{level},2), 1:a:size(gp{level},1), cos(-cls{level}.GaborOrientation(1:a:end, 1:a:end)), sin(-cls{level}.GaborOrientation(1:a:end, 1:a:end)), 'r');
    
    
%     strong_gradient = gradient_mag > (0.1*max(gradient_mag(:)));
%     aligned_gradient = abs(normal_to_gradient - cls{level}.GaborOrientation) < pi/6;
%     
%     figure; imagesc(normal_to_gradient); axis image; colormap(hsv(256));
%     figure; imagesc(cls{level}.GaborOrientation); axis image; colormap(hsv(256));
%     figure; imagesc(gradient_mag); axis image; colormap(gray(256));
%     
%     figure; imagesc(strong_gradient); axis image;
%     figure; imagesc(aligned_gradient); axis image;
%     g_mag = cls{level}.GaborResponse;
%     g_mag(strong_gradient | aligned_gradient) = 0;
%     figure; imagesc(g_mag); axis image; colormap(gray(256));
end
%%
for level = 4:5
    figure;
    imagesc(-cls{level}.GaborOrientation); axis image; colormap(hsv(256)); hold on;
    a= 1;
    quiver(1:a:size(gp{level},2), 1:a:size(gp{level},1), cos(-cls{level}.GaborOrientation(1:a:end, 1:a:end)), sin(-cls{level}.GaborOrientation(1:a:end, 1:a:end)));
end
    
%%
%Now implemented gradient throwaway in MB_CLS_SELECTION - shall we see if
%it works?


%%
Tau = 4;
Len = 8;
NumAngles = 12;

sx = Tau / (2*sqrt(2*log(2)));
sy = Len*sx;
scale = 2*ceil(sy);

xx = repmat((-scale:scale)', 1 , 2*scale+1);
yy = repmat(-scale:scale, 2*scale+1, 1);

gabor_zero = exp(-0.5*((xx.^2 / sx^2) + (yy.^2 / sy^2))).*cos(2*pi*xx / Tau)...
    / (2*pi*sx*sy);

gabor_zero = gabor_zero ./ sum(gabor_zero(:));

clear sx sy scale xx yy

angles = -90:(180/NumAngles):90 - (180/NumAngles);
gabor_theta = cell(NumAngles, 1);

figure; hold on;
for ii = 1:NumAngles

    %make orienatated Gabor filter
    theta = angles(ii);
    gabor_theta{ii} = imrotate(gabor_zero, theta, 'bilinear', 'crop');
    gabor_theta{ii} = gabor_theta{ii} ./ sum(gabor_theta{ii}(:));
    mesh(gabor_theta{ii}, 'FaceColor', 'none');
    
end
%%
for ii = 1:NumAngles
    figure; mesh(gabor_theta{ii}, 'FaceColor', 'none'); axis([0 58 0 58 -0.125 0.26]);
end
%%
for ii = 1:NumAngles
    
    display(['Angle = ', num2str(angles(ii)),...
        ' Max = ', num2str(max(gabor_theta{ii}(:))),...
        ' Min = ', num2str(min(gabor_theta{ii}(:)))]);
    
end
%%
gabor_theta2 = cell(NumAngles, 1);

figure; hold on;
for ii = 1:NumAngles

    %make orienatated Gabor filter
    theta = angles(ii);
    gabor_theta2{ii} = imrotate(gabor_zero, theta, 'bilinear', 'crop');
    mesh(gabor_theta2{ii}, 'FaceColor', 'none');
    display(['Angle = ', num2str(angles(ii)),...
        ' Max = ', num2str(max(gabor_theta2{ii}(:))),...
        ' Min = ', num2str(min(gabor_theta2{ii}(:)))]);
    
end
%%
tic
cls_thin = cls_m.CLS > 0;

for rows = 2:511
    for cols = 2:511
        win = cls_thin(rows-1:rows+1,cols-1:cols+1);
        if sum(win(:)) > 3; cls_thin(rows, cols) = 0; end
    end
end
toc
%%
tic
cls_thin2 = cls_m.CLS > 0;
cls_counts = imfilter(double(cls_thin2), ones(3));
cls_thin2(cls_counts > 3) = 0; clear cls_counts;
toc
%%
