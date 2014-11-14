d_root = [asymmetryroot 'data/retinograms/dRIVE/test/'];
    
tp_counts = zeros(20, 102);
fp_counts = zeros(20, 102);
t_counts = zeros(20, 1);
f_counts = zeros(20, 1);
aucci = zeros(4,2);

f1 = figure('windowstyle', 'normal');
graph(f1);
set(gca,'box','on'); 
axis equal; axis([0 1 0 1]); hold all;
title('ROC curves for vessel segmentation, DRIVE database');
xlabel('FPF');
ylabel('TPF');
legend_label = cell(4,1);
colors = lines(4);
colors = colors([4 2 1 3],:);
%
for data_type = 1:4
    for ii = 1:20

        switch data_type

            case 1
                label = 'Forest DT-CWT/G2';
                %vessel_prob = load_uint8(['Z:\asym\data\retinograms\DRIVE\test\predictions\detection\dt\rf_3/' zerostr(ii,2) '_test_class.mat']);
                vessel_prob = load_uint8(['Z:\asym\data\retinograms\DRIVE\test\images\results\87410/' zerostr(ii,2) '_test_class.mat']);

            case 2
                label = 'Forest G"';
                vessel_prob = load_uint8(['C:\isbe\asymmetry_project\data\retinograms\DRIVE\test\images\results\77728/' zerostr(ii,2) '_test_class.mat']);
                %vessel_prob = load_uint8(['Z:\asym\data\retinograms\DRIVE\test\predictions\detection\g2d\rf_3/' zerostr(ii,2) '_test_class.mat']);
            case 3
                label = 'Staal';
                vessel_prob = double(rgb2gray(imread([d_root 'predictions\staal\' zerostr(ii-1, 2) '.bmp'])))/255;
            case 4
                label = 'Niemeijer';
                vessel_prob = double(rgb2gray(imread([d_root 'predictions\niemeijer\' zerostr(ii, 2) '_PC_soft.bmp'])))/255;
                vessel_prob(1,:) = [];

        end

        v_mask = u_load([d_root 'vessel_masks\' zerostr(ii,2) '_test_v_mask.mat']);
        f_mask = u_load([d_root 'foveal_masks\' zerostr(ii,2) '_test_f_mask.mat']);

        %Compute ROC counts for image
        [~, ~, tp_count fp_count] = calculate_roc_image(vessel_prob, v_mask,(-1:100)/100, f_mask, 'dilate', 0);
        t_counts(ii) = sum(v_mask(f_mask));
        f_counts(ii) = sum(~v_mask(f_mask));

        %Increment total counts
        tp_counts(ii,:) = tp_count;
        fp_counts(ii,:) = fp_count;


    end

    aucci(data_type,:) = bootci(2000,@compute_auc, fp_counts, tp_counts, f_counts, t_counts);
        
    %Compute ROC points for complete set of data
    roc_pts = [sum(fp_counts)' / sum(f_counts) sum(tp_counts)' / sum(t_counts)];

    %Compute AUC for ROC curve
    auc = sum( (roc_pts(2:end,1)-roc_pts(1:end-1,1)) .* roc_pts(1:end-1,2)) + ...
            0.5*sum( (roc_pts(2:end,1)-roc_pts(1:end-1,1)) .* (roc_pts(2:end,2)-roc_pts(1:end-1,2)) );
        
    figure(f1);
    plot(roc_pts(:,1), roc_pts(:,2), '-', 'color', colors(data_type,:)); axis([0 1 0 1]);
    legend_label{data_type,1} = [ label ' AUC: ' num2str(auc)];%label;%
end
legend(legend_label, 'location', 'southeast');
%exportfig('Z:\asym\data\retinograms\2012 MICCAI\centre\segmentation_roc.png', 'png');
%exportfig('Z:\asym\data\retinograms\2012 MICCAI\centre\segmentation_roc.pdf');
%%
d_root = [asymmetryroot 'data/retinograms/DRIVE/test/'];
vessel_ori = load_uint8([d_root 'predictions\orientation\dt\rf_3\02_ori.mat']);
vessel_ori_g = load_uint8([d_root 'predictions\orientation\\g2d\analytic\orientations\002_test_ori.mat']);
vessel_prob = load_uint8([d_root 'predictions/detection/dt/rf_3/02_test_class.mat']);
v_mask = u_load([d_root 'vessel_masks/02_test_v_mask.mat']);
ret = u_load([d_root 'images/02_test.mat']);
%%
x_lim = [15 100];
y_lim = [275 365];

f1 = figure(...
    'windowstyle', 'normal',...
    'units', 'pixels',...
    'position', [100 100 5*(diff(x_lim)+1) 5*(diff(y_lim)+1)]);
imgray(v_mask);
a1 = gca;
set(a1,...
    'units', 'pixels',...
    'position', [0 0 5*(diff(x_lim)+1) 5*(diff(y_lim)+1)],...
    'xlim', x_lim,...
    'ylim', y_lim,...
    'xtick', [], 'xticklabel', [],...
    'ytick', [], 'yticklabel', []);

f2 = figure(...
    'windowstyle', 'normal',...
    'units', 'pixels',...
    'position', [100 100 5*(diff(x_lim)+1) 5*(diff(y_lim)+1)]);
imgray(vessel_prob);
a2 = gca;
set(a2,...
    'units', 'pixels',...
    'position', [0 0 5*(diff(x_lim)+1) 5*(diff(y_lim)+1)],...
    'xlim', x_lim,...
    'ylim', y_lim,...
    'xtick', [], 'xticklabel', [],...
    'ytick', [], 'yticklabel', []);

f3 = figure(...
    'windowstyle', 'normal',...
    'units', 'pixels',...
    'position', [100 100 5*(diff(x_lim)+1) 5*(diff(y_lim)+1)]);
image(ret); axis image; hold all;
a3 = gca;
set(a3,...
    'units', 'pixels',...
    'position', [0 0 5*(diff(x_lim)+1) 5*(diff(y_lim)+1)],...
    'xlim', x_lim,...
    'ylim', y_lim,...
    'xtick', [], 'xticklabel', [],...
    'ytick', [], 'yticklabel', []);

xyi = [36 297;  42 307; 57 291; 41 317; 49 344; 37 350; 68 346; 84 345; 80 325];%40 302;
theta = linspace(0, 2*pi);

for m = 1:size(xyi, 1);
    mu = angle(vessel_ori(xyi(m,2), xyi(m,1)));
    rho = abs(vessel_ori(xyi(m,2), xyi(m,1)));
    sigma = sqrt(-2*log(rho));
    display(sigma);
    r = pdf('norm', mb_mod(2*(theta - mu), 2*pi), 0, 2*sigma);
    x = 250*r.*cos(theta) / sum(r);
    y = 250*r.*sin(theta) / sum(r);
    plot(a1, xyi(m,1)+x, xyi(m,2)-y, 'linewidth', 3);
    plot(a2, xyi(m,1)+x, xyi(m,2)-y, 'linewidth', 3);
    plot(a3, xyi(m,1)+x, xyi(m,2)-y, 'linewidth', 3);
    %plot(xyi(m,1)+x, xyi(m,2)-y, '.');
end
figure(f1); exportfig('K:\isbe\conferences_and_symposia\rank_prize_fund\figures\track_oris_v_mask.pdf');
figure(f2); exportfig('K:\isbe\conferences_and_symposia\rank_prize_fund\figures\track_oris_v_prob.pdf');
figure(f3); exportfig('K:\isbe\conferences_and_symposia\rank_prize_fund\figures\track_oris_ret_rgb.pdf');
close(f1);
close(f2);
close(f3);
%%
[yc xc] = find(bwmorph(v_mask, 'thin', 'inf'));
discard = (yc < 275) | (yc > 365) | (xc < 15) | (xc > 100);
yc(discard) = [];
xc(discard) = [];

figure(...
    'windowstyle', 'normal',...
    'units', 'pixels',...
    'position', [100 100 5*(diff(x_lim)+1) 5*(diff(y_lim)+1)]);
image(ret); axis image; hold all;
set(gca,...
    'units', 'pixels',...
    'position', [0 0 5*(diff(x_lim)+1) 5*(diff(y_lim)+1)],...
    'xlim', x_lim,...
    'ylim', y_lim,...
    'xtick', [], 'xticklabel', [],...
    'ytick', [], 'yticklabel', []);
exportfig('K:\isbe\conferences_and_symposia\rank_prize_fund\figures\ret_roi.pdf');
plot(xc, yc, 'c.');
exportfig('K:\isbe\conferences_and_symposia\rank_prize_fund\figures\ret_roi_v_mask.pdf');
close(gcf);

figure(...
    'windowstyle', 'normal',...
    'units', 'pixels',...
    'position', [100 100 5*(diff(x_lim)+1) 5*(diff(y_lim)+1)]);
imgray(v_mask);
set(gca,...
    'units', 'pixels',...
    'position', [0 0 5*(diff(x_lim)+1) 5*(diff(y_lim)+1)],...
    'xlim', x_lim,...
    'ylim', y_lim,...
    'xtick', [], 'xticklabel', [],...
    'ytick', [], 'yticklabel', []);
exportfig('K:\isbe\conferences_and_symposia\rank_prize_fund\figures\v_mask_roi.pdf');
close(gcf);

figure(...
    'windowstyle', 'normal',...
    'units', 'pixels',...
    'position', [100 100 5*(diff(x_lim)+1) 5*(diff(y_lim)+1)]);
image(complex2rgb(vessel_ori.^2)); axis image;
set(gca,...
    'units', 'pixels',...
    'position', [0 0 5*(diff(x_lim)+1) 5*(diff(y_lim)+1)],...
    'xlim', x_lim,...
    'ylim', y_lim,...
    'xtick', [], 'xticklabel', [],...
    'ytick', [], 'yticklabel', []);
exportfig('K:\isbe\conferences_and_symposia\rank_prize_fund\figures\ori_roi.pdf');
close(gcf);

f1 = figure(...
    'windowstyle', 'normal',...
    'units', 'pixels',...
    'position', [100 100 5*(diff(x_lim)+1) 5*(diff(y_lim)+1)]);
imgray(v_mask); 
a1 = gca; 
set(a1,...
    'units', 'pixels',...
    'position', [0 0 5*(diff(x_lim)+1) 5*(diff(y_lim)+1)],...
    'xlim', x_lim,...
    'ylim', y_lim,...
    'xtick', [], 'xticklabel', [],...
    'ytick', [], 'yticklabel', []);
f2 = figure(...
    'windowstyle', 'normal',...
    'units', 'pixels',...
    'position', [100 100 5*(diff(x_lim)+1) 5*(diff(y_lim)+1)]);
image(ret); axis image; hold all; 
a2 = gca; 
set(a2,...
    'units', 'pixels',...
    'position', [0 0 5*(diff(x_lim)+1) 5*(diff(y_lim)+1)],...
    'xlim', x_lim,...
    'ylim', y_lim,...
    'xtick', [], 'xticklabel', [],...
    'ytick', [], 'yticklabel', []);
for m = 1:size(yc, 1);
    mu_rf = angle(vessel_ori(yc(m), xc(m)));
    plot(a1, xc(m)+[-4 4]*cos(mu_rf), yc(m)-[-4 4]*sin(mu_rf), 'g', 'linewidth', 2);
    plot(a2, xc(m)+[-4 4]*cos(mu_rf), yc(m)-[-4 4]*sin(mu_rf), 'g', 'linewidth', 2);
end
figure(f1); exportfig('K:\isbe\conferences_and_symposia\rank_prize_fund\figures\v_mask_roi_oris_rf.pdf');
figure(f2); exportfig('K:\isbe\conferences_and_symposia\rank_prize_fund\figures\ret_roi_oris_rf.pdf');
close(f1);
close(f2);
%%
f1 = figure(...
    'windowstyle', 'normal',...
    'units', 'pixels',...
    'position', [100 100 5*(diff(x_lim)+1) 5*(diff(y_lim)+1)]);
imgray(v_mask); 
a1 = gca; 
set(a1,...
    'units', 'pixels',...
    'position', [0 0 5*(diff(x_lim)+1) 5*(diff(y_lim)+1)],...
    'xlim', x_lim,...
    'ylim', y_lim,...
    'xtick', [], 'xticklabel', [],...
    'ytick', [], 'yticklabel', []);
f2 = figure(...
    'windowstyle', 'normal',...
    'units', 'pixels',...
    'position', [100 100 5*(diff(x_lim)+1) 5*(diff(y_lim)+1)]);
image(ret); axis image; hold all; 
a2 = gca; 
set(a2,...
    'units', 'pixels',...
    'position', [0 0 5*(diff(x_lim)+1) 5*(diff(y_lim)+1)],...
    'xlim', x_lim,...
    'ylim', y_lim,...
    'xtick', [], 'xticklabel', [],...
    'ytick', [], 'yticklabel', []);
for m = 1:size(yc, 1);
    mu_g2 = vessel_ori_g(yc(m), xc(m));
    plot(a1, xc(m)+[-4 4]*cos(mu_g2), yc(m)-[-4 4]*sin(mu_g2), 'g', 'linewidth', 2);
    plot(a2, xc(m)+[-4 4]*cos(mu_g2), yc(m)-[-4 4]*sin(mu_g2), 'g', 'linewidth', 2);
end
figure(f1); exportfig('M:\asymmetry_project\rank_prize_fund\figures\v_mask_roi_oris_g2.pdf');
figure(f2); exportfig('M:\asymmetry_project\rank_prize_fund\figures\ret_roi_oris_g2.pdf');
close(f1);
close(f2);
%%
s = load('C:\isbe\asymmetry_project\data\retinograms\DRIVE\test\mixed_orientations\02_ori_idx.mat');
xyi = [51 346;  87 346; 42 307];
theta = linspace(0, 2*pi);

figure; imgray(v_mask); a1 = gca; axis([15 100 275 365]);
for m = 1:size(xyi, 1);
    idx = s.mixed_idx(xyi(m,2), xyi(m,1));
    if idx
        mu = angle(s.mixed_oris{idx}) / 2;
        x = 5*cos(mu)';
        y = 5*sin(mu)';
        plot(a1, [xyi(m,1)-x; xyi(m,1)+x], [xyi(m,2)+y; xyi(m,2)-y], 'linewidth', 3);
    end
end
%%
s = load('C:\isbe\asymmetry_project\data\retinograms\DRIVE\test\mixed_orientations\02_ori_idx.mat');
xyi = [51 346;  87 346; 42 307];
theta = linspace(0, 2*pi);

figure; imgray(v_mask); a1 = gca; axis([15 100 275 365]);
for m = 1:size(xyi, 1);
    idx = s.mixed_idx(xyi(m,2), xyi(m,1));
    if idx
        mu = angle(s.mixed_oris{idx}) / 2;
        r = zeros(size(theta));
        for ii = 1:length(mu)
            r = r + pdf('norm', mb_mod(2*(theta - mu(ii)), 2*pi), 0, pi/3);
        end

        x = 250*r.*cos(theta) / sum(r);
        y = 250*r.*sin(theta) / sum(r);
        plot(a1, xyi(m,1)+x, xyi(m,2)-y, 'linewidth', 3);
    end
end
%%
%-------
color_map = hot(256);

spic_class = imread('K:\isbe\conferences_and_symposia\ipmi2011\spic_prob_mass046_b.bmp');
spic_im = imread('K:\isbe\conferences_and_symposia\ipmi2011\mass046.bmp');
H = rgb2hsv(spic_class);
spic_class = H(:,:,1);
spic_hot = ind2rgb(ceil(256*spic_class), color_map);

norm_class = imread('K:\isbe\conferences_and_symposia\ipmi2011\spic_prob_norm068_b.bmp');
norm_im = imread('K:\isbe\conferences_and_symposia\ipmi2011\norm068.bmp');
H = rgb2hsv(norm_class);
norm_class = H(:,:,1);
norm_hot = ind2rgb(ceil(256*norm_class), color_map);



figure; 
subplot(1,2,1); imgray(0.5*spic_hot + 0.5*double(spic_im)/256);
subplot(1,2,2); imgray(0.5*norm_hot + 0.5*double(norm_im)/256);

imwrite(spic_im, 'M:\asymmetry_project\rank_prize_fund\figures\mass046.bmp');
imwrite(norm_im, 'M:\asymmetry_project\rank_prize_fund\figures\norm068.bmp');
imwrite(0.5*spic_hot + 0.5*double(spic_im)/256, 'M:\asymmetry_project\rank_prize_fund\figures\spic_prob_mass046.bmp');
imwrite(0.5*norm_hot + 0.5*double(norm_im)/256, 'M:\asymmetry_project\rank_prize_fund\figures\spic_prob_norm068.bmp');

imwrite(spic_hot, 'M:\asymmetry_project\rank_prize_fund\figures\spic_prob_mass046_a.bmp');
imwrite(norm_hot, 'M:\asymmetry_project\rank_prize_fund\figures\spic_prob_norm068_a.bmp');

imwrite(ind2rgb(repmat((1:256)',1,10), color_map), 'M:\asymmetry_project\rank_prize_fund\figures\heat_bar.bmp');