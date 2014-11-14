warning('off', 'load_uint8:missing_variables');
mkdir C:\isbe\asymmetry_project\data\theta_masks\2004_screening_processed\mass_roi\

px_per_mm = 100/9;
sigma_range = [3.2 4.0 5.4 6.2 7.8];
r_max = 400;
r_min = px_per_mm*5;
R = px_per_mm*4;

xx = repmat(-400:399, 800, 1);
yy = repmat((399:-1:-400)', 1, 800);

theta = mod(atan2(yy, xx), pi);
rij = sqrt(xx.^2 + yy.^2);
phi = abs(asin(R./rij));

pij = 2*R./(pi*rij);

offset_r = 0;
offset_c = 0;
spacing = 2;

%load in the lists of mammograms names
abnormal_names = u_load('C:\isbe\asymmetry_project\data\mam_names\2004_screening_abnormals.mat');

for kk = 1:146;
 
    ii = kk;
    display(['processing mammo ' num2str(ii)]);
    
    %load in the f1 map and mask
    f1 = load_uint8(['C:\isbe\asymmetry_project\data\k_stellate_maps\g2d_max_scale'...
        '\2004_screening_processed\abnormals\' abnormal_names{ii} '_f1.mat']);
    mask = u_load(['C:\isbe\asymmetry_project\data\k_stellate_maps\g2d_max_scale'...
        '\2004_screening_processed\abnormals\' abnormal_names{ii} '_mask.mat']);

    %Convert f1 vector into map and smooth
    mask = mask(offset_r + (1:spacing:end), offset_c + (1:spacing:end));
    f1_map = zeros(size(mask));
    f1_map(mask) = f1;
    f1_map = imfilter(f1_map, fspecial('gaussian', 21, 4));
    
    %load in mass border
    mass_xy = u_load(['C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\meta\'...
        abnormal_names{ii} '_meta.mat']);
    [r c] = size(f1_map);
    x = mass_xy(:,1) * c;
    y = mass_xy(:,2) * r;

    %Workout location of maximum f1 score within mass border
    mass_mask = poly2mask(x,y,r,c);
    [yc xc] = find(mass_mask);
    [dummy mass_max_idx] = max(f1_map(mass_mask));
    xc = spacing*xc(mass_max_idx);
    yc = spacing*yc(mass_max_idx);
    
    clear f1_map mask; pack;

    %Load in regions associated with mass
    mam_mask = u_load(['C:\isbe\asymmetry_project\data\masks\2004_screening_processed\abnormals\'...
        abnormal_names{ii} '_mask.mat']);
    mam_mask = sample_window(mam_mask, 800, yc, xc, 0);
    
    rf_pm = load_uint8(['C:\isbe\asymmetry_project\data\orientation_maps\rf_thin\2004_screening_processed\abnormals\'...
        abnormal_names{ii} '_class.mat']);
    rf_pm = sample_window(rf_pm, 800, yc, xc, 0);

    rt_ori = mod(angle(rf_pm),pi);
    rt_line = abs(rf_pm);
    
    pack;
    
    %Look at target pixels
    rt_theta_diff = abs(theta - rt_ori);
    rt_theta_mask = ...
        ((rt_theta_diff < phi) | (abs(rt_theta_diff-pi) < phi))...
        & (rij > r_min) ...
        & (rij < r_max) ...
        & mam_mask;

    %Save the mask
    save(['C:\isbe\asymmetry_project\data\masks\2004_screening_processed\theta_roi\' ...
        abnormal_names{ii} '_mask.mat'], 'rt_theta_mask');


end
%%
warning('off', 'load_uint8:missing_variables');
offset_r = 0;
offset_c = 0;
spacing = 2;

%load in the lists of mammograms names
abnormal_names = u_load('C:\isbe\asymmetry_project\data\mam_names\2004_screening_abnormals.mat');
mkdir C:\isbe\asymmetry_project\data\mammograms\2004_screening\theta_roi
for kk = 1:146;
 
    ii = kk;
    display(['processing mammo ' num2str(ii)]);
    
    %load in the f1 map and mask
    f1 = load_uint8(['C:\isbe\asymmetry_project\data\k_stellate_maps\g2d_max_scale'...
        '\2004_screening_processed\abnormals\' abnormal_names{ii} '_f1.mat']);
    mask = u_load(['C:\isbe\asymmetry_project\data\k_stellate_maps\g2d_max_scale'...
        '\2004_screening_processed\abnormals\' abnormal_names{ii} '_mask.mat']);

    %Convert f1 vector into map and smooth
    mask = mask(offset_r + (1:spacing:end), offset_c + (1:spacing:end));
    f1_map = zeros(size(mask));
    f1_map(mask) = f1;
    f1_map = imfilter(f1_map, fspecial('gaussian', 21, 4));
    
    %load in mass border
    mass_xy = u_load(['C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\meta\'...
        abnormal_names{ii} '_meta.mat']);
    [r c] = size(f1_map);
    x = mass_xy(:,1) * c;
    y = mass_xy(:,2) * r;

    %Workout location of maximum f1 score within mass border
    mass_mask = poly2mask(x,y,r,c);
    [yc xc] = find(mass_mask);
    [dummy mass_max_idx] = max(f1_map(mass_mask));
    xc = spacing*xc(mass_max_idx);
    yc = spacing*yc(mass_max_idx);
    
    clear f1_map mask; pack;

    %Load in regions associated with mass
    roi = load_uint8(['C:\isbe\asymmetry_project\data\mammograms\2004_screening_processed\abnormals\'...
        abnormal_names{ii} '.mat']);
    mam_mask = u_load(['C:\isbe\asymmetry_project\data\masks\2004_screening_processed\abnormals\'...
        abnormal_names{ii} '_mask.mat']);
    out_val = mean(roi(~mam_mask));
    mam_mask = sample_window(mam_mask, 800, yc, xc, 0);
    
    roi = sample_window(roi, 800, yc, xc, 0);
    roi(~mam_mask) = out_val;
    %figure; imagesc(roi); axis image; colormap(gray(256));
    save(['C:\isbe\asymmetry_project\data\mammograms\2004_screening\theta_roi\'...
        abnormal_names{ii} '_roi.mat']);

end
%%
for ii = 1:146
    s = load(['C:\isbe\asymmetry_project\data\mammograms\2004_screening\theta_roi\'...
        abnormal_names{ii} '_roi.mat']);
    roi = s.roi;
    save(['C:\isbe\asymmetry_project\data\mammograms\2004_screening\theta_roi\'...
        abnormal_names{ii} '_roi.mat'], 'roi');
end
%%
for ii = 11:20
    mask = u_load(['C:\isbe\asymmetry_project\data\masks\2004_screening_processed\theta_roi\' ...
        abnormal_names{ii} '_mask.mat']);
    roi = u_load(['C:\isbe\asymmetry_project\data\mammograms\2004_screening_processed\theta_roi\'...
        abnormal_names{ii} '_roi.mat']);
    
    roi = (roi - min(roi(:))) / (max(roi(:))-min(roi(:)));
    
    roi_rgb = cat(3, 0.8*roi + 0.2*mask, 0.8*roi, 0.8*roi);
    
    figure; image(roi_rgb); axis image;
end
%%
warning('off', 'load_uint8:missing_variables');
mkdir C:\isbe\asymmetry_project\data\misc\k_locations\
offset_r = 0;
offset_c = 0;
spacing = 2;

%load in the lists of mammograms names
abnormal_names = u_load('C:\isbe\asymmetry_project\data\mam_names\2004_screening_abnormals.mat');

for ii = 1:146;
    
    display(['processing mammo ' num2str(ii)]);
    
    %load in the f1 map and mask
    f1 = load_uint8(['C:\isbe\asymmetry_project\data\k_stellate_maps\g2d_max_scale'...
        '\2004_screening_processed\abnormals\' abnormal_names{ii} '_f1.mat']);
    mask = u_load(['C:\isbe\asymmetry_project\data\k_stellate_maps\g2d_max_scale'...
        '\2004_screening_processed\abnormals\' abnormal_names{ii} '_mask.mat']);

    %Convert f1 vector into map and smooth
    mask = mask(offset_r + (1:spacing:end), offset_c + (1:spacing:end));
    f1_map = zeros(size(mask));
    f1_map(mask) = f1;
    f1_map = imfilter(f1_map, fspecial('gaussian', 21, 4));
    
    %load in mass border
    mass_xy = u_load(['C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\meta\'...
        abnormal_names{ii} '_meta.mat']);
    [r c] = size(f1_map);
    x = mass_xy(:,1) * c;
    y = mass_xy(:,2) * r;

    %Workout location of maximum f1 score within mass border
    mass_mask = poly2mask(x,y,r,c);
    [yc xc] = find(mass_mask);
    [dummy mass_max_idx] = max(f1_map(mass_mask));
    xc = spacing*xc(mass_max_idx);
    yc = spacing*yc(mass_max_idx);

    %Save the mask
    save(['C:\isbe\asymmetry_project\data\misc\k_locations\' ...
        abnormal_names{ii} '_xy.mat'], 'xc', 'yc');


end
%%
px_per_mm = 100/9;
sigma_range = [3.2 4.0 5.4 6.2 7.8];
r_max = 400;
r_min = px_per_mm*5;
R = px_per_mm*4;

xx = repmat(-400:399, 800, 1);
yy = repmat((399:-1:-400)', 1, 800);

theta = mod(atan2(yy, xx), pi);
rij = sqrt(xx.^2 + yy.^2);
phi = abs(asin(R./rij));
r_mask = rij < r_max;

pij = 2*R./(pi*rij);
abnormal_names = u_load('C:\isbe\asymmetry_project\data\mam_names\2004_screening_abnormals.mat');
aa_type = {'abnormals\', 'normals\'};
for ii = 2:10
    
    s = load(['C:\isbe\asymmetry_project\data\misc\k_locations\' ...
        abnormal_names{ii} '_xy.mat'], 'xc', 'yc');
    xc = s.xc;
    yc = s.yc;
    
    mam_name{1} = abnormal_names{ii};
    if abnormal_names{ii}(4) == 'R'
        mam_name{2} = [mam_name{1}(1:3) 'L' mam_name{1}(5:6)];
    else
        mam_name{2} = [mam_name{1}(1:3) 'R' mam_name{1}(5:6)];
    end
    
    f1 = figure;
    f2 = figure;
    for aa = 1:2

        %Load in regions associated with mass
        roi = load_uint8(['C:\isbe\asymmetry_project\data\mammograms\2004_screening_processed\' aa_type{aa}...
             mam_name{aa} '.mat']);
         
        if aa == 2
            xc = size(roi,2) - xc;
        end
        roi = sample_window(roi, 800, yc, xc, 0);
        
        mam_mask = u_load(['C:\isbe\asymmetry_project\data\masks\2004_screening_processed\' aa_type{aa}...
            mam_name{aa} '_mask.mat']);
        mam_mask = sample_window(mam_mask, 800, yc, xc, 0);

        r_sum = sum(r_mask(:) & mam_mask(:));
        
        g2_ori = load_uint8(['C:\isbe\asymmetry_project\data\orientation_maps\g2d\2004_screening_processed\' aa_type{aa}...
            mam_name{aa} '_ori.mat']);
        g2_ori = sample_window(g2_ori, 800, yc, xc, 0);

        %g2_line = load_uint8(['C:\isbe\asymmetry_project\data\line_maps\g2d\2004_screening_processed\' aa_type{aa}...
        %    mam_name{aa} '_lines.mat']);
        %g2_line = sample_window(g2_line, 800, yc, xc, 0);

        %rf_ori = mod(angle(load_uint8(['C:\isbe\asymmetry_project\data\orientation_maps\rf_prob\2004_screening_processed\' aa_type{aa}...
        %    mam_name{aa} '_class.mat'])),pi);
        %rf_ori = sample_window(rf_ori, 800, yc, xc, 0);

        %rf_line = load_uint8(['C:\isbe\asymmetry_project\data\line_maps\rf_prob\2004_screening_processed\' aa_type{aa}...
        %    mam_name{aa} '_class.mat']);
        %rf_line = sample_window(rf_line, 800, yc, xc, 0);

        rf_pm = load_uint8(['C:\isbe\asymmetry_project\data\orientation_maps\rf_thin\2004_screening_processed\' aa_type{aa}...
            mam_name{aa} '_class.mat']);
        rf_pm = sample_window(rf_pm, 800, yc, xc, 0);

        rf_ori2 = mod(angle(rf_pm),pi);
        %rf_line2 = abs(rf_pm);

        pack;

        %Look at target pixels
        g2_theta_diff = abs(theta - g2_ori);
        g2_theta_mask = ((g2_theta_diff < phi) | (abs(g2_theta_diff-pi) < phi)) & (rij > r_min) & mam_mask;

%         rf_theta_diff = abs(theta - rf_ori);
%         rf_theta_mask = ((rf_theta_diff < phi) | (abs(rf_theta_diff-pi) < phi)) & (rij > r_min) & mam_mask;

        rf_theta_diff2 = abs(theta - rf_ori2);
        rf_theta_mask2 = ((rf_theta_diff2 < phi) | (abs(rf_theta_diff2-pi) < phi)) & (rij > r_min) & mam_mask;

        roi = 0.8*(roi - min(roi(:)))/(max(roi(:))-min(roi(:)));
        roi_rgb = cat(3, 0.2*g2_theta_mask + roi, 0.2*rf_theta_mask2 + roi, roi);

        rf_cs = [...
            mean(cos(2*rf_ori2(r_mask & mam_mask))),...
            mean(sin(2*rf_ori2(r_mask & mam_mask)))];
        
        rf_mean_ori = atan2(rf_cs(2), rf_cs(1))/2;
        
        rf_vec = 400*sqrt(sum(rf_cs.^2))*[cos(rf_mean_ori) sin(rf_mean_ori)];
        
        g2_cs = [...
            mean(cos(2*g2_ori(r_mask & mam_mask))),...
            mean(sin(2*g2_ori(r_mask & mam_mask)))];
        
        g2_mean_ori = atan2(g2_cs(2), g2_cs(1))/2;
        
        g2_vec = 400*sqrt(sum(g2_cs.^2))*[cos(g2_mean_ori) sin(g2_mean_ori)];
        
        [rf_mu_hat, rf_kappa_hat] = von_mises_mle(2*rf_ori2(r_mask & mam_mask));
        [g2_mu_hat, g2_kappa_hat] = von_mises_mle(2*g2_ori(r_mask & mam_mask));
        
        [rf_sample] = mod(von_mises_sample(rf_mu_hat, rf_kappa_hat, r_sum)/2,pi);
        [g2_sample] = mod(von_mises_sample(g2_mu_hat, g2_kappa_hat, r_sum)/2,pi);
        

        figure(f1);
        subplot(1,2,aa); image(roi_rgb); axis image; hold on;
        plot(400+[0 rf_vec(1)], 400+[0 -rf_vec(2)], 'r', 'linewidth', 2);
        plot(400+[0 -rf_vec(1)], 400+[0 rf_vec(2)], 'r', 'linewidth', 2);
        plot(400+[0 g2_vec(1)], 400+[0 -g2_vec(2)], 'g', 'linewidth', 2);
        plot(400+[0 -g2_vec(1)], 400+[0 g2_vec(2)], 'g', 'linewidth', 2);
        
        figure(f2);
        [goat rho] = rose(rf_ori2(r_mask & mam_mask), 60);
        smooth_r = medfilt1(rho(2:4:end), 5);
        rho(2:4:end) = smooth_r;
        rho(3:4:end) = smooth_r;
        subplot(2,2,aa); 
        polar(goat, rho); hold on;
        [goat rho] = rose(rf_sample, 60);
        polar(goat, rho, 'r:');
        
        [goat rho] = rose(g2_ori(r_mask & mam_mask), 60);
        %smooth_r = medfilt1(rho(2:4:720), 5);
        %rho(2:4:720) = smooth_r;
        %rho(3:4:720) = smooth_r;
        subplot(2,2,aa+2); 
        polar(goat, rho); hold on;
        [goat rho] = rose(g2_sample, 60);
        polar(goat, rho, 'r:');
%         figure(f1);
%         subplot(2,2,aa); image(roi_rgb); axis image;
%         
%         subplot(2,2,aa+2); 
%         [goat rho] = rose(2*rf_ori2(r_mask & mam_mask), 60);
%         smooth_r = medfilt1(rho(2:4:end), 5);
%         rho(2:4:end) = smooth_r;
%         rho(3:4:end) = smooth_r;
%         polar(goat, rho); hold on;
%         [goat rho] = rose(2*g2_ori(r_mask & mam_mask), 60);
%         %smooth_r = medfilt1(rho(2:4:720), 5);
%         %rho(2:4:720) = smooth_r;
%         %rho(3:4:720) = smooth_r;
%         polar(goat, rho, 'r--');
    end
    
    
end
%%
ml_shape = u_load('D:\isbe\dev\location\mean_ML_shape');
cc_shape = u_load('D:\isbe\dev\location\mean_CC_shape');

figure;
subplot(1,2,1); plot(ml_shape(:,1), ml_shape(:,2)); axis equal image;
subplot(1,2,2); plot(cc_shape(:,1), cc_shape(:,2)); axis equal image;
%%
cc_shape = u_load('D:\isbe\dev\location\mean_CC_shape');
lcc_names = get_mammo_info(dir('C:\isbe\asymmetry_project\data\orientation_maps\rf_small\2004_screening_processed\normals\*LCC*'));
mean_orientation_map = build_orientation_model(cc_shape, lcc_names, 'plot', 0);
save C:\isbe\asymmetry_project\data\misc\mean_orientation_map.mat mean_orientation_map
smooth_orientation_map = imfilter(mean_orientation_map, fspecial('gaussian', 41, 8));
figure; 
subplot(1,2,1); image(complex2rgb(mean_orientation_map)); axis image;
subplot(1,2,2); image(complex2rgb(smooth_orientation_map)); axis image;
%%
od1 = 'C:\isbe\asymmetry_project\data\orientation_maps\rf_thin\2004_screening_processed\normals\';
od2 = 'C:\isbe\asymmetry_project\data\orientation_maps\rf_small\2004_screening_processed\normals\';
o_list = dir([od1 '*.mat']);
mkdir(od2);
for ii = 53:188
    ori_map = load_uint8([od1 o_list(ii).name]);
    ori_map = exp(i*2*angle(ori_map));
    ori_map = imresize(ori_map, [1024 nan]);
    save([od2 o_list(ii).name], 'ori_map');
    clear ori_map;
end
%%
%%
od1 = 'C:\isbe\asymmetry_project\data\orientation_maps\g2d\2004_screening_processed\normals\';
od2 = 'C:\isbe\asymmetry_project\data\orientation_maps\g2d_small\2004_screening_processed\normals\';
o_list = dir([od1 '*.mat']);
mkdir(od2);
for ii = 1:188
    ori_map = load_uint8([od1 o_list(ii).name]);
    ori_map = exp(i*2*ori_map);
    ori_map = imresize(ori_map, [1024 nan]);
    save([od2 o_list(ii).name], 'ori_map');
    clear ori_map;
end
%%
cc_shape = u_load('D:\isbe\dev\location\mean_CC_shape');
lcc_names = get_mammo_info(dir('C:\isbe\asymmetry_project\data\orientation_maps\rf_small\2004_screening_processed\normals\*LCC*'));
mean_orientation_map = build_orientation_model(cc_shape, lcc_names, 'plot', 0,...
    'ori_dir', 'C:\isbe\asymmetry_project\data\orientation_maps\g2d_small\2004_screening_processed\normals\');
save C:\isbe\asymmetry_project\data\misc\mean_orientation_map_g2.mat mean_orientation_map
smooth_orientation_map = imfilter(mean_orientation_map, fspecial('gaussian', 41, 8));
figure; 
subplot(1,2,1); image(complex2rgb(mean_orientation_map)); axis image;
subplot(1,2,2); image(complex2rgb(smooth_orientation_map)); axis image;
%%
cc_shape = u_load('D:\isbe\dev\location\mean_CC_shape');
cc_names = get_mammo_info(dir('C:\isbe\asymmetry_project\data\orientation_maps\rf_small\2004_screening_processed\normals\*CC*'));
mean_ori_map_g2 = build_orientation_model(cc_shape, cc_names, 'plot', 0,...
    'ori_dir', 'C:\isbe\asymmetry_project\data\orientation_maps\g2d_small\2004_screening_processed\normals\');
mean_ori_map_rf = build_orientation_model(cc_shape, cc_names, 'plot', 0,...
    'ori_dir', 'C:\isbe\asymmetry_project\data\orientation_maps\rf_small\2004_screening_processed\normals\');
save C:\isbe\asymmetry_project\data\misc\mean_orientation_map.mat mean_ori_map*
figure; 
subplot(1,2,1); image(complex2rgb(mean_ori_map_g2)); axis image;
subplot(1,2,2); image(complex2rgb(mean_ori_map_rf)); axis image;
%
smooth_ori_map_g2 = imfilter(mean_ori_map_g2, fspecial('gaussian', 41, 8));
smooth_ori_map_rf = imfilter(mean_ori_map_rf, fspecial('gaussian', 41, 8));
figure; 
subplot(1,2,1); image(complex2rgb(smooth_ori_map_g2)); axis image;
subplot(1,2,2); image(complex2rgb(smooth_ori_map_rf)); axis image;
%%
mean_ori_map_g28 = build_orientation_model(cc_shape, cc_names, 'plot', 0, 'sigma', 8,...
    'ori_dir', 'C:\isbe\asymmetry_project\data\orientation_maps\g2d_small\2004_screening_processed\normals\');
mean_ori_map_rf8 = build_orientation_model(cc_shape, cc_names, 'plot', 0, 'sigma', 8,...
    'ori_dir', 'C:\isbe\asymmetry_project\data\orientation_maps\rf_small\2004_screening_processed\normals\');
save C:\isbe\asymmetry_project\data\misc\mean_orientation_map.mat mean_ori_map*
figure; 
subplot(1,2,1); image(complex2rgb(mean_ori_map_g28)); axis image;
subplot(1,2,2); image(complex2rgb(mean_ori_map_rf8)); axis image;

smooth_ori_map_g28 = imfilter(mean_ori_map_g28, fspecial('gaussian', 41, 8));
smooth_ori_map_rf8 = imfilter(mean_ori_map_rf8, fspecial('gaussian', 41, 8));
figure; 
subplot(1,2,1); image(complex2rgb(smooth_ori_map_g28)); axis image;
subplot(1,2,2); image(complex2rgb(smooth_ori_map_rf8)); axis image;
%%
spacing = 16;
og2 = abs(smooth_ori_map_g28) .* exp(0.5*i*angle(smooth_ori_map_g28));
orf = abs(smooth_ori_map_rf8) .* exp(0.5*i*angle(smooth_ori_map_rf8));
[r c] = size(mean_ori_map_rf8);
figure; 
subplot(1,2,1); image(complex2rgb(smooth_ori_map_g28)); axis image; hold on;
quiver(1:spacing:c, 1:spacing:r,...
    real(og2(1:spacing:r, 1:spacing:c)), ...
   -imag(og2(1:spacing:r, 1:spacing:c)), 'r');
quiver(1:spacing:c, 1:spacing:r,...
   -real(og2(1:spacing:r, 1:spacing:c)), ...
    imag(og2(1:spacing:r, 1:spacing:c)), 'r');
subplot(1,2,2); image(complex2rgb(smooth_ori_map_rf8)); axis image; hold on;
quiver(1:spacing:c, 1:spacing:r,...
    real(orf(1:spacing:r, 1:spacing:c)), ...
   -imag(orf(1:spacing:r, 1:spacing:c)), 'r');
quiver(1:spacing:c, 1:spacing:r,...
   -real(orf(1:spacing:r, 1:spacing:c)), ...
    imag(orf(1:spacing:r, 1:spacing:c)), 'r');
%%
ml_shape = u_load('D:\isbe\dev\location\mean_ML_shape');
ml_names = get_mammo_info(dir('C:\isbe\asymmetry_project\data\orientation_maps\rf_small\2004_screening_processed\normals\*ML*'));
ml_ori_map_g28 = build_orientation_model(ml_shape, ml_names, 'plot', 0, 'sigma', 8, 'mlo', 1,...
    'ori_dir', 'C:\isbe\asymmetry_project\data\orientation_maps\g2d_small\2004_screening_processed\normals\');
%
ml_ori_map_rf8 = build_orientation_model(ml_shape, ml_names, 'plot', 0, 'sigma', 8, 'mlo', 1,...
    'ori_dir', 'C:\isbe\asymmetry_project\data\orientation_maps\rf_small\2004_screening_processed\normals\');
figure; 
subplot(1,2,1); image(complex2rgb(ml_ori_map_g28)); axis image;
subplot(1,2,2); image(complex2rgb(ml_ori_map_rf8)); axis image;

load C:\isbe\asymmetry_project\data\misc\mean_orientation_map.mat
save C:\isbe\asymmetry_project\data\misc\mean_orientation_map.mat mean_ori_map* ml_ori_map_*
%%
figure; 
subplot(1,2,1); imagesc(abs(ml_ori_map_g2) / length(ml_names)); axis image; colorbar;
subplot(1,2,2); imagesc(abs(ml_ori_map_rf) / length(ml_names)); axis image; colorbar;
%%
figure; 
subplot(1,2,1); imagesc(abs(mean_ori_map_g2) / length(cc_names)); axis image; colorbar;
subplot(1,2,2); imagesc(abs(mean_ori_map_rf) / length(cc_names)); axis image; colorbar;
%%
load C:\isbe\asymmetry_project\data\misc\mean_orientation_map.mat *8
smooth_ori_map_g28 = imfilter(mean_ori_map_g28, fspecial('gaussian', 41, 8));
smooth_ori_map_rf8 = imfilter(mean_ori_map_rf8, fspecial('gaussian', 41, 8));
%%
%mkdir C:\isbe\asymmetry_project\data\misc\k_locations_n;
seg_dir = 'C:\isbe\asymmetry_project\data\segmentations\2004_screening\abnormals\';
mask_dir = 'C:\isbe\asymmetry_project\data\masks\2004_screening\abnormals\';

ml_names = get_mammo_info(dir(['C:\isbe\asymmetry_project\data\mammograms\'...
    '2004_screening_processed\abnormals\*ML*']));
cc_names = get_mammo_info(dir(['C:\isbe\asymmetry_project\data\mammograms\'...
    '2004_screening_processed\abnormals\*CC*']));

ml_shape = u_load('D:\isbe\dev\location\mean_ML_shape');
cc_shape = u_load('D:\isbe\dev\location\mean_CC_shape');

for jj = 1:2
    if jj == 1
        mammo_names = ml_names;
        mean_shape = ml_shape;
    else
        mammo_names = cc_names;
        mean_shape = cc_shape;
    end
    
    for ii = 1:length(mammo_names)
        
        %load mammo mask
        mask = u_load([mask_dir mammo_names{ii} '_mask.mat']);
        [rows cols] = size(mask); 
%         figure; 
%         subplot(1,2,1); imagesc(mask); axis image; hold on; colormap gray;
        clear mask;
        
        %Load in segmentation
        seg_list = dir([seg_dir '*' mammo_names{ii} '*.mat']);
        seg = u_load([seg_dir seg_list(1).name]);
        
        %compute rescale factor
        scale_factor = seg.size(1) / rows;
        
        %load in points to find in mean shape
        load(['C:\isbe\asymmetry_project\data\misc\k_locations\' ...
            mammo_names{ii} '_xy.mat'], 'xc', 'yc');
%         plot(xc, yc, 'rx');
        
        %Scale the points
        xc = xc*scale_factor;
        yc = yc*scale_factor;

        %Swap x coordinate for right breasts
        if mammo_names{ii}(4) == 'R'
            xc = seg.size(2) - xc;    
        end
        breast_xy = [xc yc];

        %Get breast border as a standard CC/MLO shape
        if jj == 1
            [breast_shape dummy T] = get_mlo_shape(seg_list, seg_dir, 50);
            breast_shape = reshape(breast_shape, [], 2);
            breast_xy = breast_xy*T.rot - repmat(T.t,size(breast_xy,1),1);

        else
            breast_shape = reshape(get_cc_shape(seg_list, seg_dir, 50),[],2);
        end    

        %Align breast shape to the target mean shape
        [dd breast_shape_a t] = mb_procrustes(mean_shape, breast_shape);

        %Map mass border and mass centre into the space of the aligned shape
        breast_xy_a = t.b*breast_xy*t.T + repmat(t.c(1,:), size(breast_xy,1),1);

        %Now thin-plate spline warp each breast shape to the target mean shape

        %Breast shape form the src points
        s_x = breast_shape_a(:,1)';
        s_y = breast_shape_a(:,2)';

        %Mean shape forms the target displacements
        z_x = mean_shape(:,1)' - min(mean_shape(:,1));
        z_y = mean_shape(:,2)' - min(mean_shape(:,2));

        %Build geometric transform
        T = geom_alignpoints_2d([z_x; z_y], [s_x; s_y], [],...
            'transform', 'spline');    

        %Warp the breast shape, mass border and mass centre to the mean shape
        %space
        [breast_xy_w] = geom_transformpoints(breast_xy_a', T)';
        
        %Round to get new x-y coordinates in mean shape
        xcn = round(breast_xy_w(1));
        ycn = round(breast_xy_w(2));
        
        %Save the new locations
        save(['C:\isbe\asymmetry_project\data\misc\k_locations_n\' ...
            mammo_names{ii} '_xy.mat'], 'xcn', 'ycn');
        
%         subplot(1,2,2); plot(z_x, z_y); hold on; axis equal ij;
%         plot(xcn, ycn, 'rx');
    end
end
%%
seg_dir = 'C:\isbe\asymmetry_project\data\segmentations\2004_screening\abnormals\';
mask_dir = 'C:\isbe\asymmetry_project\data\masks\2004_screening\abnormals\';

ml_names = get_mammo_info(dir(['C:\isbe\asymmetry_project\data\mammograms\'...
    '2004_screening_processed\abnormals\meta\*ML*']));
cc_names = get_mammo_info(dir(['C:\isbe\asymmetry_project\data\mammograms\'...
    '2004_screening_processed\abnormals\meta\*CC*']));

ml_shape = u_load('D:\isbe\dev\location\mean_ML_shape');
cc_shape = u_load('D:\isbe\dev\location\mean_CC_shape');

if_plot = 0;

for jj = 1:2
    if jj == 1
        mammo_names = ml_names;
        mean_shape = ml_shape;
    else
        mammo_names = cc_names;
        mean_shape = cc_shape;
    end
    
    for ii = 1:length(mammo_names)
        
        source_name = mammo_names{ii};
        target_name = mammo_names{ii};
        if mammo_names{ii}(4) == 'R'
            target_name(4) = 'L';
        else
            target_name(4) = 'R';
        end 
        
        %Load mask of source to get size
        mask = u_load([mask_dir source_name '_mask.mat']);
        [rows cols] = size(mask); 
        clear mask;
        
        %Load in segmentations
        seg_source = u_load([seg_dir source_name '_segmentation.mat']);
        seg_target = u_load([seg_dir target_name '_segmentation.mat']);
        
        source_breast = seg_source.breast_border(seg_source.breast_air,:);
        target_breast = seg_target.breast_border(seg_target.breast_air,:); 
                
        %load in source points
        load(['C:\isbe\asymmetry_project\data\misc\k_locations\' ...
            source_name '_xy.mat'], 'xc', 'yc');
        
        %compute rescale factor
        scale_factor = seg_source.size(1) / rows;
        
        %Scale the source points
        xc = xc*scale_factor;
        yc = yc*scale_factor;
        
        %Compute corresponding target points
        [xyct] = select_corresponding_position(source_breast, target_breast, [xc yc], 1);
        
        if if_plot
            figure;
            subplot(1,2,1);
            plot(source_breast(:,1), source_breast(:,2)); hold on; axis equal ij;
            plot(xc, yc, 'rx');
            subplot(1,2,2);
            plot(target_breast(:,1), target_breast(:,2)); hold on; axis equal ij;
            plot(xyct(1), xyct(2), 'rx');
        end
        
        %Rescale and save the target points
        xc = round(xyct(1) / scale_factor);
        yc = round(xyct(2) / scale_factor);
        save(['C:\isbe\asymmetry_project\data\misc\k_locations\' ...
            target_name '_xy.mat'], 'xc', 'yc');
        
    end
end

%%
clear;
load C:\isbe\asymmetry_project\data\misc\mean_orientation_map.mat *8
smooth_cc_map_g2 = imfilter(mean_ori_map_g28, fspecial('gaussian', 41, 8)) / 92;
smooth_cc_map_rf = imfilter(mean_ori_map_rf8, fspecial('gaussian', 41, 8)) / 92;
smooth_ml_map_g2 = imfilter(ml_ori_map_g28, fspecial('gaussian', 41, 8)) / 96;
smooth_ml_map_rf = imfilter(ml_ori_map_rf8, fspecial('gaussian', 41, 8)) / 96;
clear *_ori_map*
pack;

von_m = 0;
px_per_mm = 100/9;
sigma_range = [3.2 4.0 5.4 6.2 7.8];
r_max = 400;
r_min = px_per_mm*5;
R = px_per_mm*4;

xx = repmat(-400:399, 800, 1);
yy = repmat((399:-1:-400)', 1, 800);

theta = mod(atan2(yy, xx), pi);
rij = sqrt(xx.^2 + yy.^2);
phi = abs(asin(R./rij));
r_mask = rij < r_max;

pij = 2*R./(pi*rij);
abnormal_names = u_load('C:\isbe\asymmetry_project\data\mam_names\2004_screening_abnormals.mat');
aa_type = {'abnormals\', 'normals\'};
for ii = 11:20
    
    mam_name{1} = abnormal_names{ii};
    if abnormal_names{ii}(4) == 'R'
        mam_name{2} = [mam_name{1}(1:3) 'L' mam_name{1}(5:6)];
    else
        mam_name{2} = [mam_name{1}(1:3) 'R' mam_name{1}(5:6)];
    end    
    
    f1 = figure;
    f2 = figure;
    for aa = 1:2
        
        
        s = load(['C:\isbe\asymmetry_project\data\misc\k_locations\' ...
            mam_name{aa} '_xy.mat'], 'xc', 'yc');
        xc = s.xc;
        yc = s.yc; clear s;

        s = load(['C:\isbe\asymmetry_project\data\misc\k_locations_n\' ...
            mam_name{aa} '_xy.mat'], 'xcn', 'ycn');
        xcn = s.xcn;
        ycn = s.ycn; clear s;

        %Load in regions associated with mass
        roi = load_uint8(['C:\isbe\asymmetry_project\data\mammograms\2004_screening_processed\' aa_type{aa}...
             mam_name{aa} '.mat']);
        roi = sample_window(roi, 800, yc, xc, 0);
        
        mam_mask = u_load(['C:\isbe\asymmetry_project\data\masks\2004_screening_processed\' aa_type{aa}...
            mam_name{aa} '_mask.mat']);
        mass_mask = sample_window(mam_mask, 800, yc, xc, 0); clear mam_mask;

        r_sum = sum(r_mask(:) & mass_mask(:));
        
        g2_ori_l = load_uint8(['C:\isbe\asymmetry_project\data\orientation_maps\g2d\2004_screening_processed\' aa_type{aa}...
            mam_name{aa} '_ori.mat']);
        g2_ori = sample_window(g2_ori_l, 800, yc, xc, 0); clear g2_ori_l;

        %g2_line = load_uint8(['C:\isbe\asymmetry_project\data\line_maps\g2d\2004_screening_processed\' aa_type{aa}...
        %    mam_name{aa} '_lines.mat']);
        %g2_line = sample_window(g2_line, 800, yc, xc, 0);

        rf_pm_l = load_uint8(['C:\isbe\asymmetry_project\data\orientation_maps\rf_thin\2004_screening_processed\' aa_type{aa}...
            mam_name{aa} '_class.mat']);
        rf_pm = sample_window(rf_pm_l, 800, yc, xc, 0); clear rf_pm_l;

        rf_ori = mod(angle(rf_pm),pi);
        %rf_line = abs(rf_pm);

        pack;

        %Look at target pixels
        g2_theta_diff = abs(theta - g2_ori);
        g2_theta_mask = ((g2_theta_diff < phi) | (abs(g2_theta_diff-pi) < phi)) & (rij > r_min) & mass_mask;

        rf_theta_diff = abs(theta - rf_ori);
        rf_theta_mask = ((rf_theta_diff < phi) | (abs(rf_theta_diff-pi) < phi)) & (rij > r_min) & mass_mask;

        roi = 0.8*(roi - min(roi(:)))/(max(roi(:))-min(roi(:)));
        roi_rgb = cat(3, 0.2*rf_theta_mask + roi, 0.2*g2_theta_mask + roi, roi);

        
        rf_cs = [...
            mean(cos(2*rf_ori(r_mask & mass_mask))),...
            mean(sin(2*rf_ori(r_mask & mass_mask)))];

        rf_mean_ori = atan2(rf_cs(2), rf_cs(1))/2;

        rf_vec = 400*sqrt(sum(rf_cs.^2))*[cos(rf_mean_ori) sin(rf_mean_ori)];

        g2_cs = [...
            mean(cos(2*g2_ori(r_mask & mass_mask))),...
            mean(sin(2*g2_ori(r_mask & mass_mask)))];

        g2_mean_ori = atan2(g2_cs(2), g2_cs(1))/2;

        g2_vec = 400*sqrt(sum(g2_cs.^2))*[cos(g2_mean_ori) sin(g2_mean_ori)];
        
        if 1 
            if mam_name{aa}(5) == 'M'
                rf_global = smooth_ml_map_rf(ycn, xcn);
                g2_global = smooth_ml_map_g2(ycn, xcn);
            else
                rf_global = smooth_cc_map_rf(ycn, xcn);
                g2_global = smooth_cc_map_g2(ycn, xcn);
            end
            rf_global = abs(rf_global) .* exp(0.5*i*angle(rf_global));
            g2_global = abs(g2_global) .* exp(0.5*i*angle(g2_global));
            rf_global = 400*[real(rf_global) imag(rf_global)];
            g2_global = 400*[real(g2_global) imag(g2_global)];
        end

        if von_m
            [rf_mu_hat, rf_kappa_hat] = von_mises_mle(2*rf_ori(r_mask & mass_mask));
            [g2_mu_hat, g2_kappa_hat] = von_mises_mle(2*g2_ori(r_mask & mass_mask));

            [rf_sample] = mod(von_mises_sample(rf_mu_hat, rf_kappa_hat, r_sum)/2,pi);
            [g2_sample] = mod(von_mises_sample(g2_mu_hat, g2_kappa_hat, r_sum)/2,pi);
        end

        figure(f1);
        subplot(1,2,aa); image(roi_rgb); axis image; hold on;
        plot(400+[0 rf_vec(1)], 400+[0 -rf_vec(2)], 'r', 'linewidth', 2);
        plot(400+[0 -rf_vec(1)], 400+[0 rf_vec(2)], 'r', 'linewidth', 2);
        plot(400+[0 g2_vec(1)], 400+[0 -g2_vec(2)], 'g', 'linewidth', 2);
        plot(400+[0 -g2_vec(1)], 400+[0 g2_vec(2)], 'g', 'linewidth', 2);
        plot(400+[0 rf_global(1)], 400+[0 -rf_global(2)], 'r--', 'linewidth', 2);
        plot(400+[0 -rf_global(1)], 400+[0 rf_global(2)], 'r--', 'linewidth', 2);
        plot(400+[0 g2_global(1)], 400+[0 -g2_global(2)], 'g--', 'linewidth', 2);
        plot(400+[0 -g2_global(1)], 400+[0 g2_global(2)], 'g--', 'linewidth', 2);
        
        figure(f2);
        [goat rho] = rose(rf_ori(r_mask & mass_mask), 60);
        smooth_r = medfilt1(rho(2:4:end), 5);
        rho(2:4:end) = smooth_r;
        rho(3:4:end) = smooth_r;
        subplot(2,2,aa); 
        polar(goat, rho); hold on;
        if von_m
            [goat rho] = rose(rf_sample, 60);
            polar(goat, rho, 'r:');
        end
        
        [goat rho] = rose(g2_ori(r_mask & mass_mask), 60);
        %smooth_r = medfilt1(rho(2:4:720), 5);
        %rho(2:4:720) = smooth_r;
        %rho(3:4:720) = smooth_r;
        subplot(2,2,aa+2); 
        polar(goat, rho); hold on;
        if von_m
            [goat rho] = rose(g2_sample, 60);
            polar(goat, rho, 'r:');
        end
%         figure(f1);
%         subplot(2,2,aa); image(roi_rgb); axis image;
%         
%         subplot(2,2,aa+2); 
%         [goat rho] = rose(2*rf_ori(r_mask & mass_mask), 60);
%         smooth_r = medfilt1(rho(2:4:end), 5);
%         rho(2:4:end) = smooth_r;
%         rho(3:4:end) = smooth_r;
%         polar(goat, rho); hold on;
%         [goat rho] = rose(2*g2_ori(r_mask & mass_mask), 60);
%         %smooth_r = medfilt1(rho(2:4:720), 5);
%         %rho(2:4:720) = smooth_r;
%         %rho(3:4:720) = smooth_r;
%         polar(goat, rho, 'r--');
    end
    
    
end
%%
mkdir C:\isbe\asymmetry_project\data\masks\2004_screening_processed\temp_roi\;
mkdir C:\isbe\asymmetry_project\data\orientation_maps\g2d\2004_screening_processed\temp_roi\;
mkdir C:\isbe\asymmetry_project\data\line_maps\g2d\2004_screening_processed\temp_roi\;
mkdir C:\isbe\asymmetry_project\data\orientation_maps\rf_thin\2004_screening_processed\temp_roi\;
mkdir C:\isbe\asymmetry_project\data\orientation_maps\rf_prob\2004_screening_processed\temp_roi\;

abnormal_names = u_load('C:\isbe\asymmetry_project\data\mam_names\2004_screening_abnormals.mat');
aa_type = {'abnormals\', 'normals\'};
%%
for ii = 1:146
    
    mam_name = cell(2,1);
    mam_name{1} = abnormal_names{ii};
    if abnormal_names{ii}(4) == 'R'
        mam_name{2} = [mam_name{1}(1:3) 'L' mam_name{1}(5:6)];
    else
        mam_name{2} = [mam_name{1}(1:3) 'R' mam_name{1}(5:6)];
    end
    
    for aa = 1:2
        
        s = load(['C:\isbe\asymmetry_project\data\misc\k_locations\' ...
            mam_name{aa} '_xy.mat'], 'xc', 'yc');
        xc = s.xc;
        yc = s.yc; clear s;

        %Load in regions associated with mass        
%         mam_mask = u_load(['C:\isbe\asymmetry_project\data\masks\2004_screening_processed\abnormals\'...
%             mam_name{aa} '_mask.mat']);
%         roi = sample_window(mam_mask, 800, yc, xc, 0); clear mam_mask; %#ok
%         save(['C:\isbe\asymmetry_project\data\masks\2004_screening_processed\temp_roi\'...
%             mam_name{aa} '_roi.mat'], 'roi');
%         
%         g2_ori = load_uint8(['C:\isbe\asymmetry_project\data\orientation_maps\g2d\2004_screening_processed\abnormals\'...
%             mam_name{aa} '_ori.mat']);
%         roi = sample_window(g2_ori, 800, yc, xc, 0); clear g2_ori; %#ok
%         save(['C:\isbe\asymmetry_project\data\orientation_maps\g2d\2004_screening_processed\temp_roi\'...
%             mam_name{aa} '_roi.mat'], 'roi');
%         
%         g2_line = load_uint8(['C:\isbe\asymmetry_project\data\line_maps\g2d\2004_screening_processed\abnormals\'...
%            mam_name{aa} '_lines.mat']);
%         roi = sample_window(g2_line, 800, yc, xc, 0); clear g2_roi; %#ok
%         save(['C:\isbe\asymmetry_project\data\line_maps\g2d\2004_screening_processed\temp_roi\'...
%             mam_name{aa} '_roi.mat'], 'roi');
% 
%         rf_pm = load_uint8(['C:\isbe\asymmetry_project\data\orientation_maps\rf_thin\2004_screening_processed\abnormals\'...
%             mam_name{aa} '_class.mat']);
%         roi = sample_window(rf_pm, 800, yc, xc, 0); clear rf_pm; %#ok
%         save(['C:\isbe\asymmetry_project\data\orientation_maps\rf_thin\2004_screening_processed\temp_roi\'...
%             mam_name{aa} '_roi.mat'], 'roi');
        
        rf_pm = load_uint8(['C:\isbe\asymmetry_project\data\orientation_maps\rf_prob\2004_screening_processed\abnormals\'...
            mam_name{aa} '_class.mat']);
        roi = sample_window(rf_pm, 800, yc, xc, 0); clear rf_pm;
        save(['C:\isbe\asymmetry_project\data\orientation_maps\rf_prob\2004_screening_processed\temp_roi\'...
            mam_name{aa} '_roi.mat'], 'roi');

        %pack;

    end
    
    
end
