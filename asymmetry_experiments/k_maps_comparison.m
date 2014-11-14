function k_maps_comparison

warning('off', 'load_uint8:missing_variables');

px_per_mm = 100/9;
sigma_range = [3.2 4.0 5.4 6.2 7.8];
r_max = 3*round(px_per_mm*sigma_range);
r_min = px_per_mm*4;
R = px_per_mm*2;

xx = repmat(-400:399, 800, 1);
yy = repmat((399:-1:-400)', 1, 800);

theta = mod(atan2(yy, xx), pi);
rij = sqrt(xx.^2 + yy.^2);
phi = abs(asin(R./rij));

pij = 2*R./(pi*rij);

offset_r = 0;
offset_c = 0;
spacing = 2;

r_sums = zeros(5,1);
r_mask = false(800,800,5);
for rr = 1:length(r_max)
    r_mask(:,:,rr) = (rij < r_max(rr)) & (rij > r_min);
    r_sums(rr) = sum(sum(r_mask(:,:,rr)));
end

%load in the lists of mammograms names
abnormal_names = u_load('C:\isbe\asymmetry_project\data\mam_names\2004_screening_abnormals.mat');

rp_u_scores = zeros(146,5);
rp_n_scores = zeros(146,5);
rp_N_scores = zeros(146,5);
rp_p_scores = zeros(146,5);

rt_u_scores = zeros(146,5);
rt_n_scores = zeros(146,5);
rt_N_scores = zeros(146,5);
rt_p_scores = zeros(146,5);

g2_u_scores = zeros(146,5);
g2_n_scores = zeros(146,5);
g2_N_scores = zeros(146,5);
g2_p_scores = zeros(146,5);

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

    g2_ori = load_uint8(['C:\isbe\asymmetry_project\data\orientation_maps\g2d\2004_screening_processed\abnormals\'...
        abnormal_names{ii} '_ori.mat']);
    g2_ori = sample_window(g2_ori, 800, yc, xc, 0);

    g2_line = load_uint8(['C:\isbe\asymmetry_project\data\line_maps\g2d\2004_screening_processed\abnormals\'...
        abnormal_names{ii} '_lines.mat']);
    g2_line = sample_window(g2_line, 800, yc, xc, 0);

    rp_ori = mod(angle(load_uint8(['C:\isbe\asymmetry_project\data\orientation_maps\rf_prob\2004_screening_processed\abnormals\'...
        abnormal_names{ii} '_class.mat'])),pi);
    rp_ori = sample_window(rp_ori, 800, yc, xc, 0);

    rp_line = load_uint8(['C:\isbe\asymmetry_project\data\line_maps\rf_prob\2004_screening_processed\abnormals\'...
        abnormal_names{ii} '_class.mat']);
    rp_line = sample_window(rp_line, 800, yc, xc, 0);
    
    rf_pm = load_uint8(['C:\isbe\asymmetry_project\data\orientation_maps\rf_thin\2004_screening_processed\abnormals\'...
        abnormal_names{ii} '_class.mat']);
    rf_pm = sample_window(rf_pm, 800, yc, xc, 0);

    rt_ori = mod(angle(rf_pm),pi);
    rt_line = abs(rf_pm);
    
    pack;
    
    %Look at target pixels
    g2_theta_diff = abs(theta - g2_ori);
    g2_theta_mask = ((g2_theta_diff < phi) | (abs(g2_theta_diff-pi) < phi)) & (rij > r_min);

    rp_theta_diff = abs(theta - rp_ori);
    rp_theta_mask = ((rp_theta_diff < phi) | (abs(rp_theta_diff-pi) < phi)) & (rij > r_min);
    
    rt_theta_diff = abs(theta - rt_ori);
    rt_theta_mask = ((rt_theta_diff < phi) | (abs(rt_theta_diff-pi) < phi)) & (rij > r_min);

    %
    rp_lm = rp_line > 0.5 & (rij > r_min);  
    rt_lm = rt_line > 0.35 & (rij > r_min);
    g2_lm = g2_line & (rij > r_min);


    for rr = 1:length(r_max)

        rp_u_scores(ii,rr) = sum(rp_theta_mask(r_mask(:,:,rr)));
        rp_n_scores(ii,rr) = sum(rp_theta_mask(r_mask(:,:,rr)) & rp_lm(r_mask(:,:,rr)));
        rp_N_scores(ii,rr) = sum(rp_lm(r_mask(:,:,rr)));
        rp_p_scores(ii,rr) = sum(rp_lm(r_mask(:,:,rr)) .* pij(r_mask(:,:,rr)));

        rt_u_scores(ii,rr) = sum(rt_theta_mask(r_mask(:,:,rr)));
        rt_n_scores(ii,rr) = sum(rt_theta_mask(r_mask(:,:,rr)) & rt_lm(r_mask(:,:,rr)));
        rt_N_scores(ii,rr) = sum(rt_lm(r_mask(:,:,rr)));
        rt_p_scores(ii,rr) = sum(rt_lm(r_mask(:,:,rr)) .* pij(r_mask(:,:,rr)));


        g2_u_scores(ii,rr) = sum(g2_theta_mask(r_mask(:,:,rr))) / r_sums(rr);
        g2_n_scores(ii,rr) = sum(g2_theta_mask(r_mask(:,:,rr)) & g2_lm(r_mask(:,:,rr)));
        g2_N_scores(ii,rr) = sum(g2_lm(r_mask(:,:,rr)));
        g2_p_scores(ii,rr) = sum(g2_lm(r_mask(:,:,rr)) .* pij(r_mask(:,:,rr)));
    end


end