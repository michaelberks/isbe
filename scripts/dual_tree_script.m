%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Don't run before you can walk... try swapping code between two dual trees
% to make sure it is feasible to modify coefficients and reconstruct an
% realistic region

region1 = imread('C:\isbe\dev\background\images\normal512\o04_001LCC_1024_3427_865.bmp');
region2 = imread('C:\isbe\dev\background\images\normal512\o04_001RCC_1024_2965_2874.bmp');

%%
dual_tree1 = dtwavexfm2(region1, 5, 'near_sym_b','qshift_b');
dual_tree2 = dtwavexfm2(region2, 5, 'near_sym_b','qshift_b');

mix_tree12 = dual_tree1;
mix_tree21 = dual_tree2;

for level = 1:4
    lim(1) = 2^(6-level) + 1;
    lim(2) = 7*2^(6-level);
    
    for ori = 1:6
        mix_tree12{level}(lim(1):lim(2),lim(1):lim(2),ori) = ...
            dual_tree2{level}(lim(1):lim(2),lim(1):lim(2),ori);
        
        mix_tree21{level}(lim(1):lim(2),lim(1):lim(2),ori) = ...
            dual_tree1{level}(lim(1):lim(2),lim(1):lim(2),ori);
    end
end

mix_image12 = dtwaveifm2(mix_tree12, 'near_sym_b','qshift_b');
mix_image21 = dtwaveifm2(mix_tree21, 'near_sym_b','qshift_b');

figure; imagesc(region1); axis image; colormap(gray(256));
figure; imagesc(region2); axis image; colormap(gray(256));
figure; imagesc(mix_image12); axis image; colormap(gray(256));
figure; imagesc(mix_image21); axis image; colormap(gray(256));
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syn_args.TargetImage = imread('C:\isbe\dev\background\images\normal512\o04_001LCC_1024_3427_865.bmp');
syn_args.FilledImage = ones(size(syn_args.TargetImage));
syn_args.FilledImage(100:199, 100:199) = 0;
syn_args.ModelDir = 'C:\isbe\dev\background\models\dual_tree\';
syn_args.ModelName = 'normal512_k10_w1_5_w2_3';
syn_args.Plot = 1;
syn_args.ConditionLevels = 1;
%%
[syn_image, dual_tree, cluster_image] = mb_gmm_dual_tree_synthesis(syn_args);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
region = ...
    imread('C:\isbe\dev\background\images\normal512\o04_001LCC_1024_3427_865.bmp');
dual_tree = dtwavexfm2(region, 5, 'near_sym_b', 'qshift_b');

syn_args.PathToTextureGMM = 'C:/isbe/dev/background/models/dual_tree/normal512_k10_w1_5_w2_3_4_3';
syn_args.SampleImage2 = dual_tree{5}(:,:,3);
syn_args.TargetImage = dual_tree{4}(:,:,3);
syn_args.FilledImage = ones(size(dual_tree{4}(:,:,3)));
syn_args.FilledImage(5:28,5:28) = 0;
syn_args.WindowSize = 5;
syn_args.WindowSize2 = 3;

[syn_image, cluster_image] = mb_gmm_tex_synthesis2_complex(syn_args);
%%
figure;
subplot(2,2,1); imagesc(abs(syn_args.TargetImage)); axis image;
subplot(2,2,2); imagesc(angle(syn_args.TargetImage)); axis image;
subplot(2,2,3); imagesc(abs(syn_image)); axis image;
subplot(2,2,4); imagesc(angle(syn_image)); axis image;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I want to look at what structure look like in the sub-bands
ii = 6;
dual_tree = u_load(['C:\isbe\dev\background\dual_tree\mass_2\mass_bg_', zerostr(ii,3), '_dual_tree.mat']);
cls = u_load(['C:\isbe\dev\background\cls\mass_2\mass_cls_', zerostr(ii,3), '.mat']);

f_handles = zeros(5,6);
for ori = 1:6
    for level = 1:5
        figure('Name', ['Level = ', num2str(level), ' Orientation = ', num2str(ori)]);
        subplot(2,2,1); imagesc(abs(dual_tree{level}(:,:,ori))); axis image;
        subplot(2,2,2); imagesc(angle(dual_tree{level}(:,:,ori))); axis image;
        subplot(2,2,3); imagesc(real(dual_tree{level}(:,:,ori))); axis image;
        subplot(2,2,4); imagesc(imag(dual_tree{level}(:,:,ori))); axis image;
    end
end
%Generate oriented CLS map
%[cls_map] =  mb_cls_map_orientations(cls);