roi2 = rot90(double(imread('C:\isbe\dev\background\images\normal512\o04_011RML_1024_3309_2896.bmp')));
load C:\isbe\dev\classification\rf\rf_reg_ori_01.mat


%roi2_prob = u_load('C:\isbe\dev\image_data\predict_masses512x512_chen\probability_image028.mat');
%roi2_prob = 1 - roi2_prob;

%
[roi2_ori] = classify_image(...
    'image_in', roi2,...
    'forest', random_forest,...
    'forest_type', 'regression',...
    'num_levels', 5);


%%
roi2_prob = rot90(1-u_load(['M:\chen\data\predict_normal_512x512\probability_image', zerostr(19,3), '.mat']));
roi2_ori = pi*roi2_ori/180;
%% 
roi2_combined = roi2_prob .* exp(i*roi2_ori);

figure; 
imagesc(roi2_ori); axis image; colormap(hsv(180));

figure; image(complex2rgb(roi2_combined.*roi2_combined)); axis image; hold on;
quiver(1:8:512, (1:8:512)', cos(roi2_ori(1:8:512,1:8:512)), -sin(roi2_ori(1:8:512,1:8:512)), 'w');

f2 = figure(...
    'windowstyle', 'normal',...
    'Units', 'pixels',...
    'position', [50 50 512 512],...
    'PaperPositionMode','auto');
axes(...
    'Units', 'pixels',...
    'position', [0 0 512 512]); 

axis ij equal off; hold on;

image(complex2rgb(roi2_combined.*roi2_combined));
print('-dtiff', '-noui', '-painters', f2, '-r300', 'K:\isbe\conferences_and_symposia\Monday Meeting 10-04-19\normal_ori.tif');
%%
