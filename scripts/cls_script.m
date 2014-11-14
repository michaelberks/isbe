%cls_script

%1. Linear detection experiment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%mkdir('C:\isbe\mammograms\new_CAD\BMP_2004_cls\');
m_list = dir('C:\isbe\mammograms\new_CAD\BMP_2004_masks\*.mat');

GaborFilterArgs.Normalise = 1;
GaborFilterArgs.HighPassFilter = true;

for ii = 11:length(m_list)
    
    %load mammogram and mask of breast region
    mammogram = double(u_load(['C:\isbe\mammograms\new_CAD\BMP_2004_half\', m_list(ii).name(1:end-9)]));
    mask = u_load(['C:\isbe\mammograms\new_CAD\BMP_2004_masks\', m_list(ii).name]);

    %build Gaussian Pyramid of mammogram
    [gp p_sizes] = buildGpyr(mammogram, 4);
    gp = mb_change_pyramid_form(gp, p_sizes, 'g'); clear p_sizes

    %Resize mask for each level of the pyramid and set non-breast region to
    %zero
    masks = cell(4,1);
    for level = 1:4
        masks{level} = imresize(~mask, size(gp{level}));
        %gp{level}(masks{level}) = 0;
    end
    clear mask;

    %Perform CLS detection on each level of the pyramid
    cls = cell(4,1);

    for level = 2:4
        cls{level} = mb_cls_selection('ImageIn', gp{level},...
            'GaborFilterArgs', GaborFilterArgs,...
            'MinLength', 12, 'Connectivity', 0, 'Thin', 1,...
            'IgnoreMap', masks{level}, 'GaborThreshold', 0,...
            'GradientThreshold', 0, 'GradientAlignment', pi/12,...
            'NMS', true, 'Debug', false);
    end
    
    %save the output
    
    save(['C:\isbe\mammograms\new_CAD\BMP_2004_cls\', m_list(ii).name(1:end-9), '_cls', 'cls']);
    display(['Finsished CLS detection in mammogram ', num2str(ii)]);
    clear masks mammogram gp cls
    
end
%%
% Let's look at the CLS detection in all the mammograms
%c_list = dir('C:\isbe\mammograms\new_CAD\BMP_2004_cls\*.mat');

for ii = 21:30
    mammogram = double(u_load(['C:\isbe\mammograms\new_CAD\BMP_2004_half\', c_list(ii).name(1:end-8)]))';
    cls = u_load(['C:\isbe\mammograms\new_CAD\BMP_2004_cls\', c_list(ii).name]);
    
    figure; imagesc(mammogram); colormap(gray(256)); axis image; hold on; caxis([10, 250]);
    
    colors = '_rgb_';
    for level = 3:4
        factor = 2^(level - 1);
        [y_pts x_pts] = find(cls{level}.CLS' > 0);
        plot(x_pts*factor, y_pts*factor, [colors(level), '.'], 'MarkerSize', 4);
    end
    
    clear mammogram cls;
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

% Do synthesised regions have less detected structures than real regions?
% Take a mammogram, select a large region in it to synthesise using the
% DT-DWT models. Compare CLS detection results from the original mammogram
% to synthesised image...

c_list = dir('C:\isbe\mammograms\new_CAD\BMP_2004_cls\*.mat');
%ii = randsample(1:length(c_list), 1);
ii =1;
%1. Load mammogram, deletcted CLS and mask, build GP of real mammo:
%mammogram_real = double(u_load(['C:\isbe\mammograms\new_CAD\BMP_2004_half\', c_list(ii).name(1:end-8)]));
mask = double(u_load(['C:\isbe\mammograms\new_CAD\BMP_2004_masks\', c_list(ii).name(1:end-8), '_mask']));
%cls_real = u_load(['C:\isbe\mammograms\new_CAD\BMP_2004_cls\', c_list(ii).name]);

%[gp_real p_sizes] = buildGpyr(double(mammogram_real), 4);
%gp_real = mb_change_pyramid_form(gp_real, p_sizes, 'g'); clear p_sizes

%%
%2. Set up arguments for synthesising new region
%fill_mask = ~roipoly(uint8(mammogram_real));
dt_args.ModelDir = [mberksroot, 'background/models/dual_tree/'];
dt_args.ModelName = 'normal512_k10_w1_3_w2_3';
dt_args.SynthesisMode = 'patch-wise';
dt_args.TargetDualTree = u_load(['C:\isbe\mammograms\new_CAD\BMP_2004_dt\', c_list(ii).name(1:end-8), '_dual_tree']);
dt_args.FilledImage = imresize(~mask, size(dt_args.TargetDualTree{1}(:,:,1)));
clear mask;
%%
%3 synthesise new region
[mammogram_syn, dual_tree, cluster_image] = mb_gmm_dual_tree_synthesis_oris2(dt_args);

%4 Build GP of synthesised mammogram and do CLS detection
[gp_syn p_sizes] = buildGpyr(double(mammogram_syn), 4);
gp_syn = mb_change_pyramid_form(gp_syn, p_sizes, 'g'); clear p_sizes

masks = cell(4,1);
for level = 1:4
    masks{level} = imresize(~mask, size(gp_syn{level}));
    %gp{level}(masks{level}) = 0;
end
clear mask;


cls_syn = cell(4,1);
%%
GaborFilterArgs.Normalise = 1;
GaborFilterArgs.HighPassFilter = true;
for level = 2:4
    cls_syn{level} = mb_cls_selection('ImageIn', gp_syn{level},...
        'GaborFilterArgs', GaborFilterArgs,...
        'MinLength', 12, 'Connectivity', 0, 'Thin', 1,...
        'IgnoreMap', masks{level}, 'GaborThreshold', 0,...
        'GradientThreshold', 0, 'GradientAlignment', pi/12,...
        'NMS', true, 'Debug', false);
end
%%

%5 Compare real and synthesised CLS
for level = 2:4
    [y_real{level} x_real{level}] = find(cls_real{level}.CLS' & ~cls_syn_full{level}.CLS');
    [y_syn{level} x_syn{level}] = find(cls_syn_full{level}.CLS' & ~cls_real{level}.CLS');
    [y_both{level} x_both{level}] = find(cls_real{level}.CLS' & cls_syn_full{level}.CLS');
    
    display(['CLS pixels in both: ', num2str(length(y_both))]);
    display(['CLS pixels in real only: ', num2str(length(y_real))]);
    display(['CLS pixels in syn only: ', num2str(length(y_syn))]);
    display([])
    
    figure; %imagesc(gp_real{level}' - gp_syn{level}'); colormap(gray(256)); axis image; hold on;
    %plot(x_real{level}, y_real{level}, 'r.', 'MarkerSize', 4);
    plot(x_syn{level}, y_syn{level}, 'b.', 'MarkerSize', 4);
    plot(x_both{level}, y_both{level}, 'b.', 'MarkerSize', 4);
end
    

     
    