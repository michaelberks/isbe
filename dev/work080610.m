dt = u_load('C:\isbe\dev\background\dual_tree\normal512\o04_001LCC_dual_tree.mat');

filled_image = ones(32);
filled_image(9:24,9:24) = 0;
syn_args.FilledImage = filled_image;
%%
syn_args.PathToTextureGMM = 'C:\isbe\dev\background\models\dual_tree\normal512_k10_w1_5_w2_1_4_1_0';
syn_args.UpperBand1 = imag(dt{4}(:,:,1));
syn_args.UpperBand2 = imag(dt{4}(:,:,2));
syn_args.UpperBand3 = imag(dt{4}(:,:,3));
syn_args.LowerBand1 = imag(dt{5}(:,:,1));
syn_args.LowerBand2 = imag(dt{5}(:,:,2));
syn_args.LowerBand3 = imag(dt{5}(:,:,3));

[imag_tree1 cluster_image] = mb_gmm_tex_synthesis_oris(syn_args);
%
syn_args.PathToTextureGMM = 'C:\isbe\dev\background\models\dual_tree\normal512_k10_w1_5_w2_1_4_4_0';
syn_args.UpperBand1 = imag(dt{4}(:,:,4));
syn_args.UpperBand2 = imag(dt{4}(:,:,5));
syn_args.UpperBand3 = imag(dt{4}(:,:,6));
syn_args.LowerBand1 = imag(dt{5}(:,:,4));
syn_args.LowerBand2 = imag(dt{5}(:,:,5));
syn_args.LowerBand3 = imag(dt{5}(:,:,6));

[imag_tree4 cluster_image] = mb_gmm_tex_synthesis_oris(syn_args);
%
syn_args.PathToTextureGMM = 'C:\isbe\dev\background\models\dual_tree\normal512_k10_w1_5_w2_1_4_1_1';
syn_args.UpperBand1 = real(dt{4}(:,:,1));
syn_args.UpperBand2 = real(dt{4}(:,:,2));
syn_args.UpperBand3 = real(dt{4}(:,:,3));
syn_args.LowerBand1 = real(dt{5}(:,:,1));
syn_args.LowerBand2 = real(dt{5}(:,:,2));
syn_args.LowerBand3 = real(dt{5}(:,:,3));

[real_tree1 cluster_image] = mb_gmm_tex_synthesis_oris(syn_args);
%
syn_args.PathToTextureGMM = 'C:\isbe\dev\background\models\dual_tree\normal512_k10_w1_5_w2_1_4_4_1';
syn_args.UpperBand1 = real(dt{4}(:,:,4));
syn_args.UpperBand2 = real(dt{4}(:,:,5));
syn_args.UpperBand3 = real(dt{4}(:,:,6));
syn_args.LowerBand1 = real(dt{5}(:,:,4));
syn_args.LowerBand2 = real(dt{5}(:,:,5));
syn_args.LowerBand3 = real(dt{5}(:,:,6));

[real_tree4 cluster_image] = mb_gmm_tex_synthesis_oris(syn_args);
%%
for ii = 1:3
    real_tree{ii} = real_tree1{ii};
    real_tree{ii+3} = real_tree4{ii};
    imag_tree{ii} = imag_tree1{ii};
    imag_tree{ii+3} = imag_tree4{ii};
end

for ori = 1:6
    figure;
    subplot(2,2,1); imagesc(real(dt{4}(:,:,ori))); colormap(jet(256)); axis image;
    subplot(2,2,2); imagesc(imag(dt{4}(:,:,ori))); colormap(jet(256)); axis image;
    subplot(2,2,3); imagesc(real_tree{ori}); colormap(jet(256)); axis image;
    subplot(2,2,4); imagesc(imag_tree{ori}); colormap(jet(256)); axis image;
end

%%
data_type = 'normal512';
start_idx = 1;
dt_list = dir([mberksroot, 'background/dual_tree/' data_type , '/*dual_tree*']);


display(['--Texture synthesis script started: ' datestr(now)]);

for ii = start_idx:start_idx+88
    
    if ii > length(dt_list); break; end;
    
    dt_args.TargetDualTree = u_load([mberksroot,...
        'background/dual_tree/', data_type, '/', dt_list(ii).name]);
    
    [rows cols] = size(dt_args.TargetDualTree{1}(:,:,1));
    row_centre = round(rows / 2);
    col_centre = round(cols / 2);

    % Make a the biggest circular mask
    m = min([rows cols]);
    rad = floor((m - 128) / 2);
    [x y] = meshgrid(1:cols, 1:rows);
    dt_args.FilledImage = (x - col_centre).^2 + (y - row_centre).^2 > rad.^2;
    clear x y m rad row_centre col_centre;

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dt_args.ModelDir = [mberksroot, 'background/models/dual_tree/'];
    dt_args.ModelName = 'normal512_k10_w1_3_w2_3';
    dt_args.SaveFile = [mberksroot, 'background/syn/dual_tree/',...
        data_type, '_', dt_args.ModelName, '_b_' zerostr(ii, 3)];

    mb_gmm_dual_tree_synthesis_oris2(dt_args);

    clear dt_args;
end
%%
%%
pyr_list = dir([mberksroot, 'background/syn/pyramid/*w2_5*']);
%%
for jj = 1:length(pyr_list)
    
    load([mberksroot, 'background/syn/pyramid/', pyr_list(jj).name]);
    [pyr p_dims] = mb_change_pyramid_form(pyramid);
    synthesised_image = mb_reconSFpyr(pyr, p_dims);
    save([mberksroot, 'background/syn/pyramid/', pyr_list(jj).name],...
        'synthesised_image', 'pyramid', 'cluster_image', 'filled_image');
    clear synthesised_image pyramid cluster_image filled_image pyr p_dims
end
%%
dt_list = dir([mberksroot, 'background/syn/dual_tree/*w1_3*']);

for jj = 61:89
    
    load([mberksroot, 'background/syn/dual_tree/', dt_list(jj).name]);
    figure('Name', ['Region ', num2str(jj)]);
    imagesc(synthesised_image); axis image; colormap(gray(256));
    clear synthesised_image pyramid cluster_image filled_image
    
end
%%
for level = 1:5
    for ori = 1:4
        figure;
        imagesc(pyramid{level+1, ori}); axis image;
    end
end
%%
for level = 1:5
    for ori = 1:6
        figure;
        subplot(1,2,1); imagesc(real(dual_tree{level}(:,:,ori))); axis image;
        subplot(1,2,2); imagesc(imag(dual_tree{level}(:,:,ori))); axis image;
    end
end
%%
pyr_list = dir([mberksroot, 'background/syn/pyramid/*w2_5*']);
dt_list = dir('C:\isbe\dev\background\syn\dual_tree\w1_5_w2_1_full\*w2_1*');

for ii = 31:60
    
    phi = linspace(0, 2*pi, 200);
    p = load([mberksroot, 'background/syn/pyramid/', pyr_list(ii).name]);
    s = load([mberksroot, 'background\syn\dual_tree\w1_5_w2_1_full\', dt_list(ii).name]);
    figure;
    subplot(1,2,1); image(p.synthesised_image); axis image; colormap(gray(256));
    hold on;
    plot(192*cos(phi)+256, 192*sin(phi)+256, 'r:');
    subplot(1,2,2); image(s.synthesised_image); axis image; colormap(gray(256));
    hold on;
    plot(128*cos(phi)+256, 128*sin(phi)+256, 'r:');
end

%%
load C:/isbe/dev/background/models/dual_tree/normal512_k10_W1_3_w2_3_1_imag
cluster_model.NumClusters = 1;
save C:/isbe/dev/background/models/dual_tree/normal512_k10_W1_3_w2_3_1_imag cluster_model
load C:/isbe/dev/background/models/dual_tree/normal512_k10_W1_3_w2_3_2_imag
cluster_model.NumClusters = 1;
save C:/isbe/dev/background/models/dual_tree/normal512_k10_W1_3_w2_3_2_imag cluster_model
load C:/isbe/dev/background/models/dual_tree/normal512_k10_W1_3_w2_3_3_imag
cluster_model.NumClusters = 1;
save C:/isbe/dev/background/models/dual_tree/normal512_k10_W1_3_w2_3_3_imag cluster_model
load C:/isbe/dev/background/models/dual_tree/normal512_k10_W1_3_w2_3_4_imag
cluster_model.NumClusters = 1;
save C:/isbe/dev/background/models/dual_tree/normal512_k10_W1_3_w2_3_4_imag cluster_model
load C:/isbe/dev/background/models/dual_tree/normal512_k10_W1_3_w2_3_1_real
cluster_model.NumClusters = 1;
save C:/isbe/dev/background/models/dual_tree/normal512_k10_W1_3_w2_3_1_real cluster_model
load C:/isbe/dev/background/models/dual_tree/normal512_k10_W1_3_w2_3_2_real
cluster_model.NumClusters = 1;
save C:/isbe/dev/background/models/dual_tree/normal512_k10_W1_3_w2_3_2_real cluster_model
load C:/isbe/dev/background/models/dual_tree/normal512_k10_W1_3_w2_3_3_real
cluster_model.NumClusters = 1;
save C:/isbe/dev/background/models/dual_tree/normal512_k10_W1_3_w2_3_3_real cluster_model
load C:/isbe/dev/background/models/dual_tree/normal512_k10_W1_3_w2_3_4_real
cluster_model.NumClusters = 1;
save C:/isbe/dev/background/models/dual_tree/normal512_k10_W1_3_w2_3_4_real cluster_model