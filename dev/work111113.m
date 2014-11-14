x = repmat((1:256)-128.5, 256, 1);
circle_mask = (x.^2 + x'.^2) < 116^2; clear x;
sigma = 2;
width = round(3*sigma);
sigmasq = sigma^2;

x	= (-width:width);
g	= exp(-0.5* (x.*x)/sigmasq);
dg	= -x/sigmasq .* g;
ddg	= (-1/sigmasq * g) - (x/sigmasq .* dg);

Ixx_means = zeros(12,1);
Ixx_sds = zeros(12,1);
Iyy_means = zeros(12,1);
Iyy_sds = zeros(12,1);
Ixy_means = zeros(12,1);
Ixy_sds = zeros(12,1);
tan_theta_means = zeros(12,1);
tan_theta_sds = zeros(12,1);
theta_means = zeros(12,1);
theta_sds = zeros(12,1);

%%
for jj = 1:36
    ori = (jj-1)*5;
    stripey = zeros(256);
    stripey_label = false(256);

    phi = pi * ori / 180;
    R = [cos(phi) -sin(phi); sin(phi) cos(phi)];
    cxy = [zeros(15,1) (16:16:240)' - 128.5]*R + 128.5;

    for ii = 1:15
        %[gbar, label] = create_gauss_bar(2, 1, ori, 256, 256, cxy(ii,1), cxy(ii,2));
        [gbar, dummy, label] = create_sin_bar(2, 1, ori, 256, 256, 0.5, cxy(ii,1), cxy(ii,2));
        stripey = stripey + gbar;
        stripey_label = stripey_label | label;
    end
    stripey_label = stripey_label & circle_mask;
    stripey_n = imnoise(imnoise(stripey, 'gaussian', 0, 0.5), 'speckle');
    %stripey_n = stripey;
%     figure;
%     subplot(1,2,1); imagesc(stripey_label); axis image; colormap(gray(256));
%     subplot(1,2,2); imagesc(stripey_n); axis image; colormap(gray(256));

    Ixx = conv2(g',ddg,stripey_n,'same'); % = Ixx
    Ixy = conv2(dg',dg,stripey_n,'same'); % = Ixy
    Iyy = conv2(ddg',g,stripey_n,'same'); % = Iyy

    Ixx = Ixx(stripey_label);
    Iyy = Iyy(stripey_label);
    Ixy = -Ixy(stripey_label);
    tan_theta = (2*Ixy) ./ (Ixx - Iyy);
    theta = .5*atan(tan_theta);
    
    cc = cos(theta).^2;
    ss = sin(theta).^2;
    s2 = sin(2*theta);
    wo1 = abs(Ixx.*cc + Iyy.*ss + Ixy.*s2);
    wo2 = abs(Ixx.*ss + Iyy.*cc - Ixy.*s2);
    swap = wo2 < wo1;
    theta(swap) = theta(swap) + pi/2;
    theta = mod(theta, pi);
    theta_errs = mb_mod(theta - phi, pi);
    
    Ixx_means(jj) = mean(Ixx);
    Ixx_sds(jj) = std(Ixx);
    Iyy_means(jj) = mean(Iyy);
    Iyy_sds(jj) = std(Iyy);
    Ixy_means(jj) = mean(Ixy);
    Ixy_sds(jj) = std(Ixy);
    tan_theta_means(jj) = mean(tan_theta);
    tan_theta_sds(jj) = std(tan_theta);
    theta_means(jj) = mean(theta);
    theta_sds(jj) = std(theta);
    
    %derivs = (Ixx - Iyy).*sin(2*theta) - 2*Ixy.*cos(2*theta);
    
    figure;
%     subplot(2,2,1); hist(Ixx, 20);
%     subplot(2,2,2); hist(Ixy, 20);
%     subplot(2,2,3); hist(tan_theta, 20);
%     %subplot(2,2,4); hist(theta, linspace(0, pi, 180));
%     subplot(2,2,4); hist(theta_errs, linspace(-pi/2, pi/2, 180)); set(gca, 'ylim', [0 1700]);
    centres = linspace(-pi/2, pi/2, 150);
    counts = hist(theta_errs, centres);
    bar(centres, counts / sum(counts)); set(gca, 'ylim', [0 .1]);
    title(['\phi = ' num2str(ori) '\circ']);
%     hist(derivs, 50);
    
    
end
%%
centres = linspace(-90, 90, 150);
for jj = 1:36
    ori_l = (jj-1)*5;
    ori_u = (jj+1)*5;
    
    idx = (ori_l <= true_oris) & (true_oris < ori_u);
    num_pts = sum(idx);
    counts3 = hist(ori_errors3(idx), centres);
    counts5 = hist(ori_errors5(idx), centres);
    figure; 
    subplot(2,1,1); bar(centres, counts3 / num_pts); set(gca, 'ylim', [0 .05]);
    subplot(2,1,2); bar(centres, counts5 / num_pts); set(gca, 'ylim', [0 .05]);
    title(['\phi = ' num2str(jj*5) '\circ']);
end
%%
g3 = zeros(1000,2);
g5 = zeros(1000,2);
s3 = zeros(1000,2);
s5 = zeros(1000,2);

contrast_range = [4 8];
width_range = [4 16];
decay_rate = 4;

for jj = 1:1000
    ori = rand * 180;
    width = rand*4 + 2;
  
    mu = (contrast_range(2) - contrast_range(1)) / (2*log(decay_rate));
    contrast = contrast_range(1) + exp_rand_sample(mu);
    
    squash = rand;
    
    phi = pi * ori / 180;
    gbar = create_gauss_bar(width/2, contrast, ori, 64, 64, 32, 32);
    sbar = create_sin_bar(width/2, contrast, ori, 64, 64, squash, 32, 32);
    
    gbar = imnoise(imnoise(gbar, 'gaussian', 0, 0.5), 'speckle', 0.02);
    sbar = imnoise(imnoise(sbar, 'gaussian', 0, 0.5), 'speckle', 0.02);
    

    [g_strength3 g_ori3] = gaussian_clover_line_exp(gbar, 2, 6);
    [g_strength5 g_ori5] = gaussian_clover_line_exp(gbar, 2, 10);
    [s_strength3 s_ori3] = gaussian_clover_line_exp(sbar, 2, 6);
    [s_strength5 s_ori5] = gaussian_clover_line_exp(sbar, 2, 10);

    g3(jj,1) = g_strength3(32,32);
    g5(jj,1) = g_strength5(32,32);
    s3(jj,1) = s_strength3(32,32);
    s5(jj,1) = s_strength5(32,32);
    g3(jj,2) = mb_mod(mod(g_ori3(32,32),pi) - phi, pi);
    g5(jj,2) = mb_mod(mod(g_ori5(32,32),pi) - phi, pi);
    s3(jj,2) = mb_mod(mod(s_ori3(32,32),pi) - phi, pi);
    s5(jj,2) = mb_mod(mod(s_ori5(32,32),pi) - phi, pi);
end
%%
for s = 4%[1 2 4 8]
    for w = 5%3:4
        clover_dir = ['Z:\asym\data\retinograms\DRIVE\test\orientations\clover_w' num2str(w) '_s' num2str(s)];
        mkdir([clover_dir '\orientations\theta_a\']);
        mkdir([clover_dir '\orientations\theta_b\']);
        mkdir([clover_dir '\lines\response_a\']);
        mkdir([clover_dir '\lines\response_b\']);
        for ii = 1:20
            %load(['C:\isbe\asymmetry_project\data\synthetic_lines\real512\image' zerostr(ii,3) '.mat']);
            ret = u_load(['Z:\asym\data\retinograms\DRIVE\test\images_extended\' zerostr(ii,2) '_test_ext.mat']);
            ret = rgb2gray(ret);
            [line_map ori_map] = gaussian_clover_line_exp(ret, s, w*s);
            save([clover_dir '\orientations\' zerostr(ii,3) '_class.mat'], 'ori_map');
            save([clover_dir '\lines\' zerostr(ii,3) '_class.mat'], 'line_map');
            
        end

    end
end
%%
clover_dir = 'Z:\asym\data\retinograms\DRIVE\test\orientations\clover';
mkdir([clover_dir '\orientations\']);
mkdir([clover_dir '\lines\']);

for ii = 1:20
    %load(['C:\isbe\asymmetry_project\data\synthetic_lines\real512\image' zerostr(ii,3) '.mat']);
    ret = u_load(['Z:\asym\data\retinograms\DRIVE\test\images_extended\' zerostr(ii,2) '_test_ext.mat']);
    ret = rgb2gray(ret);
    [line_map ori_map] = gaussian_clover_line(ret, [1 2 4 8]);
    save([clover_dir '\orientations\' zerostr(ii,3) '_ori.mat'], 'ori_map');
    save([clover_dir '\lines\' zerostr(ii,3) '_line.mat'], 'line_map');

end
%%
for s = 4%[1 2 4 8]
    for w = 3:5
        clover_dir = ['Z:\asym\data\retinograms\DRIVE\test\orientations\clover_w' num2str(w) '_s' num2str(s)];
        
        v_oris_gt = [];
        v_oris_pred = [];
        for ii = 1:20
            edge_mask = logical(imread(['Z:\asym\data\retinograms\DRIVE\test\mask\' zerostr(ii,2) '_test_mask.gif']));
            vessel_mask = logical(imread(['Z:\asym\data\retinograms\DRIVE\test\1st_manual\' zerostr(ii,2) '_manual1.gif']));
            vessel_mask(~edge_mask) = false;
            ori_gt = u_load(['Z:\asym\data\retinograms\DRIVE\test\orientations\' zerostr(ii,2) '_ori1.mat']);
            ori_map = u_load([clover_dir '\orientations\' zerostr(ii,3) '_class.mat']);
            
            v_oris_gt = [v_oris_gt; ori_gt(vessel_mask)];%#ok
            v_oris_pred = [v_oris_pred; ori_map(vessel_mask)];%#ok
            
        end
        v_oris_gt = mod(90*angle(v_oris_gt) / pi, 180);
        v_oris_pred = mod(180*v_oris_pred / pi, 180);
        v_oris_diff = mb_mod(v_oris_gt - v_oris_pred, 180);
        save([clover_dir '\v_ori_errs.mat'], 'v_oris_*');
        figure; plot(v_oris_gt, v_oris_diff, 'r.');

    end
end
%%
v3 = load('Z:\asym\data\retinograms\DRIVE\test\orientations\clover_w3_s4\v_ori_errs.mat');
v5 = load('Z:\asym\data\retinograms\DRIVE\test\orientations\clover_w5_s4\v_ori_errs.mat');
v_err_dist3 = hist3([v3.v_oris_gt v3.v_oris_diff], [45 180]);
v_err_dist5 = hist3([v5.v_oris_gt v5.v_oris_diff], [45 180]);

figure;
subplot(2,1,1); imagesc(v_err_dist3');
subplot(2,1,2); imagesc(v_err_dist5');
%%  
clover_dir = 'Z:\asym\data\retinograms\DRIVE\test\images_extended\results\313846\';
v_ori_c_gt = [];
v_ori_c_pred = [];
ori_c_errs = [];
v_ori_gt = [];
v_ori_pred = [];
ori_errs = [];
for ii = 1:20
    edge_mask = logical(imread(['Z:\asym\data\retinograms\DRIVE\test\mask\' zerostr(ii,2) '_test_mask.gif']));
    vessel_mask = u_load(['Z:\asym\data\retinograms\DRIVE\test\vessel_masks\' zerostr(ii,2) '_test_v_mask.mat']);
    vessel_mask(~edge_mask) = false;
    
    vessel_c_mask = bwmorph(vessel_mask, 'skel', 'inf');
    
    ori_gt = u_load(['Z:\asym\data\retinograms\DRIVE\test\orientations\' zerostr(ii,2) '_ori1.mat']);
    ori_map = load_uint8([clover_dir zerostr(ii,2) '_test_ext_class.mat']);%_ori
    
    %ori_errs = [ori_errs; ori_error(ori_gt(vessel_mask), exp(complex(0, 2*ori_map(vessel_mask))))]; %#ok
    ori_errs = [ori_errs; ori_error(ori_gt(vessel_mask & ~vessel_c_mask), exp(complex(0, 2*angle(ori_map(vessel_mask & ~vessel_c_mask)))))]; %#ok
    
    
    %ori_c_errs = [ori_errs; ori_error(ori_gt(vessel_c_mask), exp(complex(0, 2*ori_map(vessel_c_mask))))]; %#ok
    ori_c_errs = [ori_c_errs; ori_error(ori_gt(vessel_c_mask), exp(complex(0, 2*angle(ori_map(vessel_c_mask)))))]; %#ok
    
    
%     v_ori_gt = [v_ori_gt; ori_gt(vessel_mask & ~vessel_c_mask)];%#ok
%     v_ori_pred = [v_ori_pred; ori_map(vessel_mask & ~vessel_c_mask)];%#ok
%     v_ori_c_gt = [v_ori_c_gt; ori_gt(vessel_c_mask)];%#ok
%     v_ori_c_pred = [v_ori_c_pred; ori_map(vessel_c_mask)];%#ok
    

end
% v_ori_c_gt = mod(90*angle(v_ori_c_gt) / pi, 180);
% v_ori_c_pred = mod(180*anglv_ori_c_pred / pi, 180);
% v_ori_c_diff = mb_mod(v_ori_c_gt - v_ori_c_pred, 180);
% 
% v_ori_gt = mod(90*angle(v_ori_gt) / pi, 180);
% v_ori_pred = mod(180*v_ori_pred / pi, 180);
% v_ori_diff = mb_mod(v_ori_gt - v_ori_pred, 180);
% 
% save([clover_dir 'v_ori_errs.mat'], 'v_ori_*');
% figure; plot(v_ori_c_gt, v_ori_c_diff, 'r.');
%%
%%
sigma = 2;
width = round(5*sigma);
sigmasq = sigma^2;

x	= (-width:width);
g	= exp(-0.5* (x.*x)/sigmasq);
dg	= -x/sigmasq .* g;
ddg	= (-1/sigmasq * g) - (x/sigmasq .* dg);

for ori = 1:180
    
    gbar = create_gauss_bar(2, 1, ori, 32, 32, 16, 16);
        

    Ixx = conv2(g',ddg,gbar,'same'); % = Ixx
    Ixy = conv2(dg',dg,gbar,'same'); % = Ixy
    Iyy = conv2(ddg',g,gbar,'same'); % = Iyy

    Ixx = Ixx(16,16);
    Iyy = Iyy(16,16);
    Ixy = -Ixy(16,16);
    tan_theta = (2*Ixy) ./ (Ixx - Iyy);
    theta = .5*atan(tan_theta);
    
    cc = cos(theta).^2;
    ss = sin(theta).^2;
    s2 = sin(2*theta);
    w_a = abs(Ixx.*cc + Iyy.*ss + Ixy.*s2);
    w_b = abs(Ixx.*ss + Iyy.*cc - Ixy.*s2);
    
    Ixxs(ori) = Ixx;
    Iyys(ori) = Iyy;
    Ixys(ori) = Ixy;
    w_as(ori) = w_a;
    w_bs(ori) = w_b;
   
end
%
figure;
hold on;
plot(1:180,Ixxs, 'g');
plot(1:180,Ixys, 'r');
plot(1:180,Iyys, 'b');
%%
warning('off', 'load_uint8:missing_variables');
load C:\isbe\asymmetry_project\data\misc\true_oris.mat true_oris
true_oris = pi*true_oris/180;
ori_errors3 = ...
     u_load('C:\isbe\asymmetry_project\data\synthetic_lines\real512\results\clover_w3_s4\orientations\centre_lineerrors\ori_errors.mat');
ori_errors5 = ...
     u_load('C:\isbe\asymmetry_project\data\synthetic_lines\real512\results\clover_w5_s4\orientations\centre_lineerrors\ori_errors.mat');
 
[response_a3] =...
    get_image_classifications(...
    'C:\isbe\asymmetry_project\data\synthetic_lines\real512\', ...
    'C:\isbe\asymmetry_project\data\synthetic_lines\real512\results\clover_w3_s4\lines\response_a', 'centre_line');
[response_b3] =...
    get_image_classifications(...
    'C:\isbe\asymmetry_project\data\synthetic_lines\real512\', ...
    'C:\isbe\asymmetry_project\data\synthetic_lines\real512\results\clover_w3_s4\lines\response_b', 'centre_line');

[response_a5] =...
    get_image_classifications(...
    'C:\isbe\asymmetry_project\data\synthetic_lines\real512\', ...
    'C:\isbe\asymmetry_project\data\synthetic_lines\real512\results\clover_w5_s4\lines\response_a', 'centre_line');
[response_b5] =...
    get_image_classifications(...
    'C:\isbe\asymmetry_project\data\synthetic_lines\real512\', ...
    'C:\isbe\asymmetry_project\data\synthetic_lines\real512\results\clover_w5_s4\lines\response_b', 'centre_line');

[theta_a3] =...
    get_image_classifications(...
    'C:\isbe\asymmetry_project\data\synthetic_lines\real512\', ...
    'C:\isbe\asymmetry_project\data\synthetic_lines\real512\results\clover_w3_s4\orientations\theta_a', 'centre_line');
[theta_b3] =...
    get_image_classifications(...
    'C:\isbe\asymmetry_project\data\synthetic_lines\real512\', ...
    'C:\isbe\asymmetry_project\data\synthetic_lines\real512\results\clover_w3_s4\orientations\theta_b', 'centre_line');

[theta_a5] =...
    get_image_classifications(...
    'C:\isbe\asymmetry_project\data\synthetic_lines\real512\', ...
    'C:\isbe\asymmetry_project\data\synthetic_lines\real512\results\clover_w5_s4\orientations\theta_a', 'centre_line');
[theta_b5] =...
    get_image_classifications(...
    'C:\isbe\asymmetry_project\data\synthetic_lines\real512\', ...
    'C:\isbe\asymmetry_project\data\synthetic_lines\real512\results\clover_w5_s4\orientations\theta_b', 'centre_line');

err_a3 = 180*ori_error(theta_a3, true_oris)/pi;
err_b3 = 180*ori_error(theta_b3, true_oris)/pi;
err_a5 = 180*ori_error(theta_a5, true_oris)/pi;
err_b5 = 180*ori_error(theta_b5, true_oris)/pi;
%%
idx_a = abs(ori_errors5) > 80 & abs(ori_errors3) < 80 & abs(err_b5) > 80;

figure;
subplot(2,2,1); hold on;
plot(response_a5(idx_a), response_b5(idx_a), 'r.'); plot([0 0], get(gca, 'ylim'), 'k'); plot(get(gca, 'xlim'), [0 0], 'k');
subplot(2,2,2); hold on;
plot(response_a5(~idx_a), response_b5(~idx_a), 'r.'); plot([0 0], get(gca, 'ylim'), 'k'); plot(get(gca, 'xlim'), [0 0], 'k');
subplot(2,2,3); hold on;
plot(response_a3(idx_a), response_b3(idx_a), 'b.'); plot([0 0], get(gca, 'ylim'), 'k'); plot(get(gca, 'xlim'), [0 0], 'k');
subplot(2,2,4); hold on;
plot(response_a3(~idx_a), response_b3(~idx_a), 'b.'); plot([0 0], get(gca, 'ylim'), 'k'); plot(get(gca, 'xlim'), [0 0], 'k');

figure;
subplot(2,2,1); hold on;
plot(abs(response_a5(idx_a)), abs(response_b5(idx_a)), 'r.'); plot([0 0], get(gca, 'ylim'), 'k'); plot(get(gca, 'xlim'), [0 0], 'k');
subplot(2,2,2); hold on;
plot(abs(response_a5(~idx_a)), abs(response_b5(~idx_a)), 'r.'); plot([0 0], get(gca, 'ylim'), 'k'); plot(get(gca, 'xlim'), [0 0], 'k');
subplot(2,2,3); hold on;
plot(abs(response_a3(idx_a)), abs(response_b3(idx_a)), 'b.'); plot([0 0], get(gca, 'ylim'), 'k'); plot(get(gca, 'xlim'), [0 0], 'k');
subplot(2,2,4); hold on;
plot(abs(response_a3(~idx_a)), abs(response_b3(~idx_a)), 'b.'); plot([0 0], get(gca, 'ylim'), 'k'); plot(get(gca, 'xlim'), [0 0], 'k');

figure; hist(ori_errors3(idx_a), 180);
%%
idx_b = abs(ori_errors5) > 80 & abs(ori_errors3) < 80 & abs(err_a5) > 80;

figure;
subplot(2,2,1); hold on;
plot(response_a5(idx_b), response_b5(idx_b), 'r.'); plot([0 0], get(gca, 'ylim'), 'k'); plot(get(gca, 'xlim'), [0 0], 'k');
subplot(2,2,2); hold on;
plot(response_a5(~idx_b), response_b5(~idx_b), 'r.'); plot([0 0], get(gca, 'ylim'), 'k'); plot(get(gca, 'xlim'), [0 0], 'k');
subplot(2,2,3); hold on;
plot(response_a3(idx_b), response_b3(idx_b), 'b.'); plot([0 0], get(gca, 'ylim'), 'k'); plot(get(gca, 'xlim'), [0 0], 'k');
subplot(2,2,4); hold on;
plot(response_a3(~idx_b), response_b3(~idx_b), 'b.'); plot([0 0], get(gca, 'ylim'), 'k'); plot(get(gca, 'xlim'), [0 0], 'k');

figure;
subplot(2,2,1); hold on;
plot(abs(response_a5(idx_b)), abs(response_b5(idx_b)), 'r.'); plot([0 0], get(gca, 'ylim'), 'k'); plot(get(gca, 'xlim'), [0 0], 'k');
subplot(2,2,2); hold on;
plot(abs(response_a5(~idx_b)), abs(response_b5(~idx_b)), 'r.'); plot([0 0], get(gca, 'ylim'), 'k'); plot(get(gca, 'xlim'), [0 0], 'k');
subplot(2,2,3); hold on;
plot(abs(response_a3(idx_b)), abs(response_b3(idx_b)), 'b.'); plot([0 0], get(gca, 'ylim'), 'k'); plot(get(gca, 'xlim'), [0 0], 'k');
subplot(2,2,4); hold on;
plot(abs(response_a3(~idx_b)), abs(response_b3(~idx_b)), 'b.'); plot([0 0], get(gca, 'ylim'), 'k'); plot(get(gca, 'xlim'), [0 0], 'k');

figure; hist(ori_errors3(idx_b), 180);
%%
figure;
subplot(2,2,1); hold on; axis equal;
plot(response_a5(idx_b), response_a3(idx_b), 'r.'); plot([0 0], get(gca, 'ylim'), 'k'); plot(get(gca, 'xlim'), [0 0], 'k');
xlabel('R_a, 5\sigma'); ylabel('R_a, 3\sigma');
subplot(2,2,2); hold on; axis equal;
plot(response_a5(idx_b), response_b3(idx_b), 'r.'); plot([0 0], get(gca, 'ylim'), 'k'); plot(get(gca, 'xlim'), [0 0], 'k');
xlabel('R_a, 5\sigma'); ylabel('R_b, 3\sigma');
subplot(2,2,3); hold on; axis equal;
plot(response_b5(idx_b), response_b3(idx_b), 'b.'); plot([0 0], get(gca, 'ylim'), 'k'); plot(get(gca, 'xlim'), [0 0], 'k');
xlabel('R_b, 5\sigma'); ylabel('R_b, 3\sigma');
subplot(2,2,4); hold on; axis equal;
plot(response_b5(~idx_b), response_a3(~idx_b), 'b.'); plot([0 0], get(gca, 'ylim'), 'k'); plot(get(gca, 'xlim'), [0 0], 'k');
xlabel('R_b, 5\sigma'); ylabel('R_a, 3\sigma');
%%
figure;
subplot(2,2,1); hold on;
plot(abs(response_a5(idx_b)), abs(response_b5(idx_b)), 'r.'); plot([0 0], get(gca, 'ylim'), 'k'); plot(get(gca, 'xlim'), [0 0], 'k');
subplot(2,2,2); hold on;
plot(abs(response_a5(~idx_b)), abs(response_b5(~idx_b)), 'r.'); plot([0 0], get(gca, 'ylim'), 'k'); plot(get(gca, 'xlim'), [0 0], 'k');
subplot(2,2,3); hold on;
plot(abs(response_a3(idx_b)), abs(response_b3(idx_b)), 'b.'); plot([0 0], get(gca, 'ylim'), 'k'); plot(get(gca, 'xlim'), [0 0], 'k');
subplot(2,2,4); hold on;
plot(abs(response_a3(~idx_b)), abs(response_b3(~idx_b)), 'b.'); plot([0 0], get(gca, 'ylim'), 'k'); plot(get(gca, 'xlim'), [0 0], 'k');
%%
ori_errors3(:,2) = max(abs(response_a3), abs(response_b3));
ori_errors5(:,2) = max(abs(response_a5), abs(response_b5));

errors3 = sortrows(ori_errors3, -2);
errors5 = sortrows(ori_errors5, -2);

ori_errors3s = sort(abs(ori_errors3));
ori_errors5s = sort(abs(ori_errors5));
num_pts = size(ori_errors3, 1);
%%
ori_cdf3 = zeros(101,1);
mean_pct3 = zeros(100,1);
ori_cdf5 = zeros(101,1);
mean_pct5 = zeros(100,1);
for jj = 1:100
    x = ceil(jj*num_pts/100);
    ori_cdf3(jj+1) = ori_errors3s(x,1);
    mean_pct3(jj) = mean(abs(errors3(1:x,1)));
    ori_cdf5(jj+1) = ori_errors5s(x,1);
    mean_pct5(jj) = mean(abs(errors5(1:x,1)));
end

figure;
subplot(1,2,1); hold on;
plot(ori_cdf3, linspace(0,1,101), 'r', 'linewidth', 2);
plot(ori_cdf5, linspace(0,1,101), 'g', 'linewidth', 2);
subplot(1,2,2); hold on;
plot(linspace(0,1,100), mean_pct3, 'r', 'linewidth', 2);
plot(linspace(0,1,100), mean_pct5, 'g', 'linewidth', 2);
%%
mkdir Z:\asym\data\retinograms\DRIVE\test\vessel_masks\centre_idx
for ii = 1:20
    foveal_mask = u_load(['Z:\asym\data\retinograms\DRIVE\test\foveal_masks\' zerostr(ii,2) '_test_f_mask']);
    vessel_mask = u_load(['Z:\asym\data\retinograms\DRIVE\test\vessel_masks\' zerostr(ii,2) '_test_v_mask.mat']);
    vessel_mask(~foveal_mask) = false;
    
    vessel_c_mask = bwmorph(vessel_mask, 'skel', 'inf');
    
    vessel_centres = vessel_c_mask(vessel_mask);
    save(...
        ['Z:\asym\data\retinograms\DRIVE\test\vessel_masks\centre_idx\' zerostr(ii,2) '_test_centre_idx.mat'],...
        'vessel_centres');
    

end
%%
emap = [1 1 1; 0 0 1; 0 1 0; 1 0 0; 1 1 0];
retroot = 'Z:\asym\data\retinograms\DRIVE\test\predictions\';
for ii = 16%1:20
    
    response_map_rf = load_uint8([retroot 'dt\rf_3\'  zerostr(ii,2) '_ori.mat']);
    response_map_g2 = load_uint8([retroot 'g2d\analytic\responses\'  zerostr(ii,3) '_test_response.mat']);
    error_map_rf = 180*u_load([retroot 'dt\rf_3\errors\'  zerostr(ii,2) '_error_map.mat']) / pi;
    error_map_g2 = 180*u_load([retroot 'g2d\analytic\errors\'  zerostr(ii,2) '_error_map.mat']) / pi;
    
    %Generate HSV representation for this part
    c_hsv1 = ones(size(error_map_g2,1), size(error_map_g2,2), 3);
    c_hsv1(:,:,1) = (90 + abs(error_map_g2) - abs(error_map_rf)) / 360;
    c_hsv1(:,:,3) = 1 - abs(error_map_rf / 90);
    
    c_hsv2 = ones(size(error_map_g2,1), size(error_map_g2,2), 3);
    c_hsv2(:,:,1) = (90 + abs(error_map_rf) - abs(error_map_g2)) / 360;
    c_hsv2(:,:,3) = 1 - abs(error_map_g2 / 90);

    %Convert to rgb
    c_rgb1 = hsv2rgb(c_hsv1);
    c_rgb2 = hsv2rgb(c_hsv2);
    
    %error_diff = 180*(abs(error_map_g2) - abs(error_map_rf))/pi;
    error_diffq = zeros(size(error_map_g2));
    error_diffq(~isnan(error_map_rf)) = 2;
    error_diffq((abs(error_map_rf) > 25) & (abs(error_map_g2) < 15)) = 1;
    error_diffq((abs(error_map_g2) > 25) & (abs(error_map_rf) < 15)) = 3;
    error_diffq((abs(error_map_g2) < 15) & (abs(error_map_rf) < 15)) = 4;
    
    %figure; imagesc(error_diff); axis image;% colorbar;
    %figure; imagesc(abs(response_map_rf)); axis image; a1 = gca;
    %figure; imagesc(response_map_g2); axis image; a2 = gca;
    figure; imagesc(error_diffq); axis image; colormap(emap); a3 = gca;
    figure; 
    subplot(1,2,1); image(c_rgb1); axis image;
    subplot(1,2,2); image(c_rgb2); axis image;
    %figure; imagesc(180*abs(error_map_rf)/pi); caxis([0 45]); axis image; a4 = gca; colorbar;
    
    %imwrite(uint8(error_diffq), emap, [retroot 'misc_experiments\dt_g2_comparison\' zerostr(ii,2) '_error_diff.png']); 
    
    %linkaxes([a1 a2 a3 a4]);
end
%%
ii = 16;
ret = rgb2gray(load_uint8(['Z:\asym\data\retinograms\DRIVE\test\images_extended\' zerostr(ii,2) '_test_ext']));
ret_roi = ret(100:350, 300:end);
rf_ori = angle(load_uint8([retroot 'dt\rf_3\'  zerostr(ii,2) '_ori.mat']));
rf_ori = rf_ori(100:350, 300:end);
ret_mask = u_load(['Z:\asym\data\retinograms\DRIVE\test\vessel_masks\' zerostr(ii,2) '_test_v_mask']);
ret_mask = ret_mask(100:350, 300:end);
fov_mask = u_load(['Z:\asym\data\retinograms\DRIVE\test\foveal_masks\' zerostr(ii,2) '_test_f_mask']);
fov_mask = fov_mask(100:350, 300:end);
%
%compute 1st deriv.
[g1d_r g1d_o] = gaussian_1st_derivative_gradient2(ret_roi, 2);%[1 2 4 8]

%compute 1st deriv.
[g2d_r g2d_o] = gaussian_clover_line(ret_roi, 2);%[1 2 4 8]

ridge = mb_curvature_sign_change(g1d_r, g2d_r, g2d_o, 0, 0);
ridge(~fov_mask) = 0;

ridge_rf = mb_curvature_sign_change(g1d_r, g2d_r, rf_ori, 0, 0);
ridge_rf(~fov_mask) = 0;
%
figure; 
subplot(1,2,1); imagesc(ret_roi); axis image; colormap(gray(256));
subplot(1,2,2); imagesc(ridge > 0); axis image;

figure; 
subplot(1,2,1); imagesc(ridge > 0); axis image; colormap(gray(256));
subplot(1,2,2); imagesc(ridge_rf > 0); axis image;
%
ridge_skel = bwmorph(ridge > 0, 'skel', inf);
[groups_map_g2 group_vecs_g2] = contour_linking(g2d_o, ridge_skel);
[groups_map_rf group_vecs_rf] = contour_linking(rf_ori, ridge_skel);

% figure; 
% subplot(1,2,1); imagesc(groups_map_g2); axis image;
% subplot(1,2,2); imagesc(groups_map_rf); axis image;
%%
ridge_skel = bwmorph(ret_mask, 'skel', inf);
ridge_skel(rand(size(ridge_skel)) > .5) = 0;
figure; imagesc(ridge_skel); axis image; colormap(gray(256));

[groups_map_g2 group_vecs_g2] = contour_linking(g2d_o, ridge_skel);
[groups_map_rf group_vecs_rf] = contour_linking(rf_ori, ridge_skel);

figure; 
subplot(1,2,1); imagesc(ridge_skel); colormap(gray(256)); axis image; hold on;
for jj = 1:length(group_vecs_g2)
    plot(group_vecs_g2{jj}(:,1), group_vecs_g2{jj}(:,2), 'r');
end

subplot(1,2,2); imagesc(ridge_skel); colormap(gray(256)); axis image; hold on;
for jj = 1:length(group_vecs_rf)
    plot(group_vecs_rf{jj}(:,1), group_vecs_rf{jj}(:,2), 'r');
end

figure; 
subplot(1,2,1); imagesc(ridge_skel); colormap(gray(256)); axis image; hold on;
for jj = 1:length(group_vecs_g2)
    plot(group_vecs_g2{jj}(:,1), group_vecs_g2{jj}(:,2), 'r');
end

subplot(1,2,2); imagesc(ridge_skel); colormap(gray(256)); axis image; hold on;
for jj = 1:length(group_vecs_rf)
    plot(group_vecs_rf{jj}(:,1), group_vecs_rf{jj}(:,2), 'r');
end

figure; 
subplot(1,2,1); imagesc(false(size(ret_mask))); colormap(gray(256)); axis image; hold on;
for jj = 1:length(group_vecs_g2)
    plot(group_vecs_g2{jj}(:,1), group_vecs_g2{jj}(:,2), 'r');
end

subplot(1,2,2); imagesc(false(size(ret_mask))); colormap(gray(256)); axis image; hold on;
for jj = 1:length(group_vecs_rf)
    plot(group_vecs_rf{jj}(:,1), group_vecs_rf{jj}(:,2), 'r');
end
%%
[orientation_errors mean_image_errors line_contrasts_centre] =...
    compute_image_orientation_errors(...
    'C:\isbe\asymmetry_project\data\synthetic_lines\real512', ...
    'C:\isbe\asymmetry_project\data\synthetic_lines\real512\results\44412',...
    'centre_line');
[orientation_errors mean_image_errors line_contrasts_all] =...
    compute_image_orientation_errors(...
    'C:\isbe\asymmetry_project\data\synthetic_lines\real512', ...
    'C:\isbe\asymmetry_project\data\synthetic_lines\real512\results\44412',...
    'all_line');
%%
load Z:\asym\data\retinograms\DRIVE\test\retinogram_properties.mat
load Z:\asym\data\retinograms\DRIVE\test\predictions\dt\rf_3\errors\orientation_errors.mat
vessel_widths = [];
vessel_centres = [];
for ii = 1:20
    vessel_widths = [vessel_widths; line_widths{ii}];
    vessel_centres = logical([vessel_centres; centre_inds{ii}]);
end
%%
centre_errs = prediction_errs(vessel_centres);
centre_widths = vessel_widths(vessel_centres);

centre_errs = abs(centre_errs)*180/pi;


thin_vessels = centre_widths <= 2;
thin_errs = centre_errs(thin_vessels);


num_pts = length(centre_errs);
ori_cdf = zeros(101,1);
sorted_ori_errors = sort(centre_errs);
for jj = 1:100
    x = ceil(jj*num_pts/100);
    ori_cdf(jj+1) = sorted_ori_errors(x,1);
end    
figure; plot(ori_cdf, (0:100)/100, 'linewidth', 2);

num_pts = length(thin_errs);
ori_cdf = zeros(101,1);
sorted_ori_errors = sort(thin_errs);
for jj = 1:100
    x = ceil(jj*num_pts/100);
    ori_cdf(jj+1) = sorted_ori_errors(x,1);
end    
figure; plot(ori_cdf, (0:100)/100, 'linewidth', 2);
%%
retroot = [asymmetryroot('shared'),'data\synthetic_lines\real512\'];
load([retroot,'orientations\all_gt_orientations.mat']);
load([retroot, 'contrasts\centre_gt_contrasts.mat']);

pred_type{1}  = 'analytic';             pred_decomp{1}  = 'mono';   label{1}  = 'Monogenic';
pred_type{2}  = 'analytic';             pred_decomp{2}  = 'g1d';    label{2}  = 'Gauss. 1^{st} deriv';
pred_type{3}  = 'analytic';             pred_decomp{3}  = 'g2d';    label{3}  = 'Gauss. 2^{nd} deriv';
pred_type{4}  = 'rf_3';                 pred_decomp{4}  = 'dt';     label{4}  = '3x3 Forest DT-CWT';
pred_type{5}  = 'rf_3';                 pred_decomp{5}  = 'mono';   label{5}  = '3x3 Forest Mono.';
pred_type{6}  = 'rf_3';                 pred_decomp{6}  = 'g1d';    label{6}  = '3x3 Forest G''';
pred_type{7}  = 'rf_3';                 pred_decomp{7}  = 'g2d';    label{7}  = '3x3 Forest G"';
pred_type{8}  = 'boosted_regression_3';   pred_decomp{8}  = 'dt';     label{8}  = 'DT-CWT/Boost Reg, 3x3';
pred_type{9}  = 'boosted_regression_3';   pred_decomp{9}  = 'mono';   label{9}  = 'Monogenic/Boost Reg, 3x3';
pred_type{10} = 'boosted_regression_3';   pred_decomp{10} = 'g1d';    label{10} = 'Gauss. 1^{st}/Boost Reg, 3x3';
pred_type{11} = 'boosted_regression_3';   pred_decomp{11} = 'g2d';    label{11} = 'Gauss. 2^{nd}/Roost Reg, 3x3';
pred_type{12} = 'linear_regression_3';    pred_decomp{12} = 'dt';     label{12} = 'DT-CWT/Lin Reg, 3x3';
pred_type{13} = 'linear_regression_3';    pred_decomp{13} = 'mono';   label{13} = 'Monogenic/Lin Reg, 3x3';
pred_type{14} = 'linear_regression_3';    pred_decomp{14} = 'g1d';    label{14} = 'Gauss. 1^{st}/Lin Reg, 3x3';
pred_type{15} = 'linear_regression_3';    pred_decomp{15} = 'g2d';    label{15} = 'Gauss. 2^{nd}/Lin Reg, 3x3';

% plot_num = 1;
% contrasts_x = [];
% contrasts_y = [];
% for ii = [4 7 12 3]
%     
%     err_dir = [retroot 'predictions\' pred_decomp{ii} '\' pred_type{ii} '\errors\'];
%     
%     %save the predictions and errors
%     load([err_dir 'orientation_errors.mat'], 'centre_errs');
%     centre_errs = 180*abs(centre_errs)/pi;
%     
%     
%     [contrasts_x(plot_num,:), contrasts_y(plot_num,:)] = ...
%         smoother(line_contrasts_centre, centre_errs, 100);
%         
%     plot_num = plot_num + 1;
% 
%     
%     
% end
% 
% graph(1); clf; hold on;
% plot(contrasts_x',contrasts_y');
% axis([0,max(line_contrasts_centre)*1.05,0,100]);
% set(gca,'box','on'); 
% xlabel('Line contrast'); ylabel('Abs. angular error (degrees)');
% legend(label([4 7 12 3]));
% exportfig([retroot 'figs\contrast_vs_error_mammo.pdf']);

plot_num = 1;
responses_x = [];
responses_y = [];
for ii = [3 6 7 5 4]
    
    err_dir = [retroot 'predictions\' pred_decomp{ii} '\' pred_type{ii} '\errors\'];
    
    %save the predictions and errors
    %load([err_dir 'orientation_predictions.mat'], 'predicted_orientations');
    load([err_dir 'orientation_responses.mat'], 'orientation_responses');
    load([err_dir 'orientation_errors.mat'], 'prediction_errs');
    prediction_errs = 180*abs(prediction_errs)/pi;

    if ii == 3
        orientation_responses = abs(orientation_responses / max(orientation_responses));
    end
    
    [responses_x(plot_num,:), responses_y(plot_num,:)] = ...
        kernel_smoother(orientation_responses, prediction_errs, 100);
        
    plot_num = plot_num + 1;
  
end

graph(2); clf; hold on;
plot(responses_x',responses_y');
axis([0,1,0,45]);
set(gca,'box','on'); 
xlabel('Magnitude of predicted orientation vector'); ylabel('Abs. angular error (degrees)');
legend(label([3 6 7 5 4]));
%exportfig([retroot 'figs\response_vs_error_mammo.pdf']);










































    