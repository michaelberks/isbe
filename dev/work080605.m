% Work 5th June 2008

% Test using mb_compute_decomp_corr across a group of masses

%work 29/05/2008
pyr_list = dir('C:\isbe\dev\background\pyramid\mass_2\*pyr*');

p_sum_pyr = zeros(20);
map_sum_pyr = zeros(20);

for jj = 1:length(pyr_list)
    
    display(['processing region ', num2str(jj)]);
    
    pyr = u_load(['C:\isbe\dev\background\pyramid\mass_2\', pyr_list(jj).name]);
    
    [decomp_corr] = mb_compute_decomp_corr(...
        pyr, 'pyr', 'ComputeDecompLocal', 0);

    p_sum_pyr = p_sum_pyr + decomp_corr.subbands.p_val;
    map_sum_pyr = map_sum_pyr + decomp_corr.subbands.map;
end
p_sum_pyr = p_sum_pyr / jj;
%%
figure; imagesc(map_sum_pyr); axis image;
set(gca, 'xtick', [0:4:20], 'ytick', [0:4:20]);
title('Correlation map for real steerable pyramid data')

% saveas(gcf,
% 'K:\isbe\reports\background_texture\figures\corr_map_sum_pyr.eps', 'eps');
%%
write_im_from_colormap(map_sum_pyr,...
    'K:\isbe\reports\background_texture\figures\corr_map_sum_pyr.bmp', ...
    jet(256));
%%
pyr_list = dir('C:\isbe\dev\background\syn\pyramid\normal512*');

p_sum_syn_pyr = zeros(20);
map_sum_syn_pyr = zeros(20);

for jj = 1:length(pyr_list)
    
    display(['processing region ', num2str(jj)]);
    
    syn = load(['C:\isbe\dev\background\syn\pyramid\', pyr_list(jj).name]);
    
    [decomp_corr] = mb_compute_decomp_corr(...
        syn.pyramid, 'pyr', 'ComputeDecompLocal', 0, 'SampleMask', ~syn.filled_image);

    p_sum_syn_pyr = p_sum_syn_pyr + decomp_corr.subbands.p_val;
    map_sum_syn_pyr = map_sum_syn_pyr + decomp_corr.subbands.map;
end
p_sum_syn_pyr = p_sum_syn_pyr / jj;

save C:\isbe\dev\background\structures\corr_pyr_w1_11_w2_5_subband *sum_syn_pyr

%%
figure; imagesc(map_sum_syn_pyr); axis image;
set(gca, 'xtick', [0:6:60], 'ytick', [0:6:60]);
title('Correlation map for synthetic steerable pyramid data')

write_im_from_colormap(map_sum_syn_pyr,...
    'K:\isbe\reports\background_texture\figures\corr_map_sum_syn_pyr.bmp', ...
    jet(256));
%%
dt_list = dir('C:\isbe\dev\background\dual_tree\mass_2\*dual_tree*');

p_sum_dt = zeros(60);
map_sum_dt = zeros(60);

for jj = 1:length(dt_list)
    
    display(['processing region ', num2str(jj)]);
    
    dual_tree = u_load(['C:\isbe\dev\background\dual_tree\mass_2\', dt_list(jj).name]);
    
    [decomp_corr] = mb_compute_decomp_corr(...
        dual_tree, 'dual_tree', 'ComputeDecompLocal', 0);

    p_sum_dt = p_sum_dt + decomp_corr.subbands.p_val;
    map_sum_dt = map_sum_dt + decomp_corr.subbands.map;
end
p_sum_dt = p_sum_dt / jj;
%%
figure; imagesc(map_sum_dt); axis image;
set(gca, 'xtick', [0:6:60], 'ytick', [0:6:60]);
title('Correlation map for real dual-tree data')
% saveas(gcf,
% 'K:\isbe\reports\background_texture\figures\corr_map_sum_dt.eps', 'eps');
%%
write_im_from_colormap(map_sum_dt,...
    'K:\isbe\reports\background_texture\figures\corr_map_sum_dt.bmp', ...
    jet(256));
%%


dt_list = dir('C:\isbe\dev\background\syn\dual_tree\*w1_3*');

p_sum_syn_dt = zeros(60);
map_sum_syn_dt = zeros(60);

for jj = 1:length(dt_list)
    
    display(['processing region ', num2str(jj)]);
    
    syn = load(['C:\isbe\dev\background\syn\dual_tree\', dt_list(jj).name]);
    
    [decomp_corr] = mb_compute_decomp_corr(...
        syn.dual_tree, 'dual_tree', 'ComputeDecompLocal', 0, 'SampleMask', ~syn.filled_image);
    
    clear syn;
    
    p_sum_syn_dt = p_sum_syn_dt + decomp_corr.subbands.p_val;
    map_sum_syn_dt = map_sum_syn_dt + decomp_corr.subbands.map;
    
    clear decomp_corr;
end
p_sum_syn_dt = p_sum_syn_dt / jj;


save C:\isbe\dev\background\structures\corr_dt_w1_3_w2_3_full_subband *sum_syn_dt

%%
figure; imagesc(map_sum_syn_dt); axis image;
set(gca, 'xtick', [0:6:60], 'ytick', [0:6:60]);
title('Correlation map for synthetic dual-tree data (model: W1 = 5, W2 = 0)')

%%
% saveas(gca,
% 'K:\isbe\reports\background_texture\figures\corr_map_syn_sum_dt_0.eps');
write_im_from_colormap(map_sum_syn_dt,...
    'K:\isbe\reports\background_texture\figures\corr_map_syn_sum_dt_0.bmp', ...
    jet(256));
%%

pyr_list = dir('C:\isbe\dev\background\pyramid\normal512\*pyr*');
cls_list = dir('C:\isbe\dev\background\cls\normal512\*cls*');

p_sum_pyr_cls = zeros(20);
map_sum_pyr_cls = zeros(20);
p_sum_pyr_non = zeros(20);
map_sum_pyr_non = zeros(20);


for jj = 1:length(pyr_list)
    
    display(['processing region ', num2str(jj)]);
    
    cls = u_load(['C:\isbe\dev\background\cls\normal512\', cls_list(jj).name]);
    pyr = u_load(['C:\isbe\dev\background\pyramid\normal512\', pyr_list(jj).name]);
    
    filled_image_cls = cls{1}.CLS;
    filled_image_cls(1:32, :) = 0;
    filled_image_cls(end-31:end, :) = 0;
    filled_image_cls(:,1:32) = 0;
    filled_image_cls(:,end-31:end) = 0;
    
    filled_image_non = ~cls{1}.CLS;
    filled_image_non(1:32, :) = 0;
    filled_image_non(end-31:end, :) = 0;
    filled_image_non(:,1:32) = 0;
    filled_image_non(:,end-31:end) = 0;
    
    [decomp_corr_cls] = mb_compute_decomp_corr(...
        'Decomp', pyr, 'ComputeDecompLocal', 0, 'SampleMask', filled_image_cls);
    
    [decomp_corr_non] = mb_compute_decomp_corr(...
        'Decomp', pyr, 'ComputeDecompLocal', 0, 'SampleMask', filled_image_non);

    p_sum_pyr_cls = p_sum_pyr_cls + decomp_corr_cls.subbands.p_val;
    map_sum_pyr_cls = map_sum_pyr_cls + decomp_corr_cls.subbands.map;
    
    p_sum_pyr_non = p_sum_pyr_non + decomp_corr_non.subbands.p_val;
    map_sum_pyr_non = map_sum_pyr_non + decomp_corr_non.subbands.map;
end
%%

dt_list = dir('C:\isbe\dev\background\dual_tree\normal512\*dual_tree*');
cls_list = dir('C:\isbe\dev\background\cls\normal512\*cls*');

p_sum_dt_cls = zeros(60);
map_sum_dt_cls = zeros(60);
p_sum_dt_non = zeros(60);
map_sum_dt_non = zeros(60);


for jj = 1:length(dt_list)
    
    display(['processing region ', num2str(jj)]);
    
    cls = u_load(['C:\isbe\dev\background\cls\normal512\', cls_list(jj).name]);
    dt = u_load(['C:\isbe\dev\background\dual_tree\normal512\', dt_list(jj).name]);
    
    filled_image_cls = cls{2}.CLS | imresize(cls{3}.CLS, [256 256]);
    filled_image_cls(1:32, :) = 0;
    filled_image_cls(end-31:end, :) = 0;
    filled_image_cls(:,1:32) = 0;
    filled_image_cls(:,end-31:end) = 0;
    
    filled_image_non = ~(cls{2}.CLS | imresize(cls{3}.CLS, [256 256]));
    filled_image_non(1:32, :) = 0;
    filled_image_non(end-31:end, :) = 0;
    filled_image_non(:,1:32) = 0;
    filled_image_non(:,end-31:end) = 0;
    
    [decomp_corr_cls] = mb_compute_decomp_corr(...
        dt, 'dt', 'ComputeDecompLocal', 0, 'SampleMask', filled_image_cls);
    
    [decomp_corr_non] = mb_compute_decomp_corr(...
        dt, 'dt', 'ComputeDecompLocal', 0, 'SampleMask', filled_image_non);

    p_sum_dt_cls = p_sum_dt_cls + decomp_corr_cls.subbands.p_val;
    map_sum_dt_cls = map_sum_dt_cls + decomp_corr_cls.subbands.map;
    
    p_sum_dt_non = p_sum_dt_non + decomp_corr_non.subbands.p_val;
    map_sum_dt_non = map_sum_dt_non + decomp_corr_non.subbands.map;
end
%%

p_sum_syn_pyr = zeros(20);
map_sum_syn_pyr = zeros(20);

for jj = 1:15
    
    pyr = u_load([mberksroot, 'background/syn/normal_2_2levels/conditioned_synthesis_circle_6_', zerostr(jj, 3)]);
    
    [rows cols] = size(pyr{2,1});
    row_centre = round(rows / 2);
    col_centre = round(cols / 2);

    % Make a the biggest circular mask
    m = min([rows cols]);
    rad = floor((m - 128) / 2);
    [x y] = meshgrid(1:cols, 1:rows);
    filled_image = (x - col_centre).^2 + (y - row_centre).^2 < rad.^2;
    
    [decomp_corr] = mb_compute_decomp_corr(...
        'Decomp', pyr, 'ComputeDecompLocal', 0, 'SampleMask', filled_image);
    
    p_sum_syn_pyr = p_sum_syn_pyr + decomp_corr.subbands.p_val;
    map_sum_syn_pyr = map_sum_syn_pyr + decomp_corr.subbands.map;
    
end
 
%%

pyr_list = dir('C:\isbe\dev\background\pyramid\normal512\*pyr*');

p_sum_local_pyr = zeros(20, 121);
map_sum_local_pyr = zeros(20, 121);

for jj = 1:1%length(pyr_list)
    
    display(['processing region ', num2str(jj)]);
    
    pyr = u_load(['C:\isbe\dev\background\pyramid\normal512\', pyr_list(jj).name]);
    
    [decomp_corr] = mb_compute_decomp_corr(...
        pyr, 'pyr', 'ComputeDecompSubbands', 0);

    p_sum_local_pyr = p_sum_local_pyr + decomp_corr.local.p_val;
    map_sum_local_pyr = map_sum_local_pyr + decomp_corr.local.map;
end
p_sum_local_pyr = p_sum_local_pyr / jj;
%%

figure;
ii = 0;
for level = 1:5
    for ori = 1:4
        
        ii = ii+1;
        subplot(5,4,ii); imagesc(reshape(map_sum_local_pyr_syn(ii,:),11,11)); axis image;
        set(gca, 'xtick', [], 'ytick', []);
        title(['Level = ', num2str(level), ', Orientation = ', num2str(ori)]);
    end
end

% saveas(gcf,
% 'K:\isbe\reports\background_texture\figures\corr_map_pyr_local.eps');
%%

dt_list = dir('C:\isbe\dev\background\dual_tree\normal512\*dual_tree*');

p_sum_local_dt = zeros(30, 121);
map_sum_local_dt = zeros(30, 121);

for jj = 1:length(dt_list)
    
    display(['processing region ', num2str(jj)]);
    
    dt = u_load(['C:\isbe\dev\background\dual_tree\normal512\', dt_list(jj).name]);
    
    [decomp_corr] = mb_compute_decomp_corr(...
        dt, 'dt', 'ComputeDecompSubbands', 0);

    p_sum_local_dt = p_sum_local_dt + decomp_corr.local.p_val;
    map_sum_local_dt = map_sum_local_dt + decomp_corr.local.map;
end
p_sum_local_dt = p_sum_local_dt / jj;
%%
figure;
ii = 0;
for level = 1:5
    for ori = 1:6
        
        ii = ii+1;
        subplot(5,6,ii); imagesc(reshape(map_sum_local_dt_syn(ii,:),11,11)); axis image; caxis([0,89]);
        set(gca, 'xtick', [], 'ytick', []);
        title(['Lev = ', num2str(level), ', Ori = ', num2str(ori)]);
    end
end

% saveas(gcf,
% 'K:\isbe\reports\background_texture\figures\corr_map_dt_local.eps');
%%
figure;
ii = 0;
for level = 1:5
    for ori = 1:4
        
        ii = ii+1;
        subplot(5,4,ii); imagesc(reshape(map_sum_local_pyr_syn(ii,:),11,11)); axis image; %caxis([0,89]);
        set(gca, 'xtick', [], 'ytick', []);
        title(['Lev = ', num2str(level), ', Ori = ', num2str(ori)]);
    end
end

% saveas(gcf, 'K:\isbe\reports\background_texture\figures\corr_map_pyr_syn_local.eps');

%%
figure;
ii = 0;
for level = 1:5
    for ori = 1:6
        
        ii = ii+1;
        subplot(5,6,ii); imagesc(reshape(map_sum_local_dt(ii,:),11,11)); axis image; caxis([0,89]);
        set(gca, 'xtick', [], 'ytick', []);
        title(['Lev = ', num2str(level), ', Ori = ', num2str(ori)]);
    end
end


% saveas(gcf, 'K:\isbe\reports\background_texture\figures\corr_map_dt_local11.eps');
    
    
    