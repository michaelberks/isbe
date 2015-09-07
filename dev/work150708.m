s = load('N:\nailfold Capillaroscopy\camera_capture\wellcome_nailfold_study\001wellcome\2015_02_27\L1_09_38_22\segments\segment_data.mat');
nailfold = imresize(s.segment_mosaic, 0.5);
%%
nailfold_mask = ~isnan(nailfold);
nailfold(~nailfold_mask) = mean(nailfold(nailfold_mask));

detect_capillaries(nailfold, 1:7, ...
    'nailfold_mask',            nailfold_mask,...
    'save_path',                'C:\isbe\nailfold\test_analysis2.mat',...
    'plot',                     1);
%%
d_list = dir('C:\isbe\nailfold\data\rsa_study\images\all_data\*d.jpg');
c_list = dir('C:\isbe\nailfold\data\rsa_study\images\*c.png');
%%
num_d = length(d_list);
reject_d = false(num_d,1);
d_names = cell(num_d,1);
for i_d = 1:num_d
    d_names{i_d} = d_list(i_d).name(1:5);
    reject_d(i_d) = length(d_list(i_d).name) ~= 10;
end
d_names(reject_d) = [];
num_d = length(d_names);

num_c = length(c_list);
reject_c = false(num_c,1);
c_names = cell(num_c,1);
for i_c = 1:num_c
    c_names{i_c} = c_list(i_c).name(1:5);
    reject_c(i_c) = length(c_list(i_c).name) ~= 10;
end
c_names(reject_c) = [];
num_c = length(c_names);
%%
shared_names = intersect(c_names, d_names);

d_im = imread(['C:\isbe\nailfold\data\rsa_study\images\all_data\' shared_names{1} 'd.jpg']);
c_im = imread(['C:\isbe\nailfold\data\rsa_study\images\' shared_names{1} 'c.png']);

%%
%derm_id_data = get_image_id_data('save_path', 'derm_image_id_data.mat'); 
derm_id_data = u_load('C:\isbe\nailfold\data\rsa_study\data_lists\derm_image_id_data.mat'); 
%%
visit = 1;
hand = 'L';
digit = 4;

c_dir = 'C:\isbe\nailfold\data\rsa_study\images\';
d_dir = 'C:\isbe\nailfold\data\rsa_study\images\all_data\';
    
for id = 1:20
    
    match = ...
        derm_id_data.people_id == id &...
        derm_id_data.visit == visit &...
        strcmpi(derm_id_data.hand, hand) &...
        derm_id_data.digit == digit;
    
    if any(match)
        
        c_name = derm_id_data.im_names{match};
        d_name = derm_id_data.derm_names{match};
        
        if exist([c_dir c_name '.png'], 'file') && ...
            ~isempty(d_name) && ...
            exist([d_dir d_name{1} 'd.jpg'], 'file')  
        
            c_im = imread([c_dir c_name '.png']);
            d_im = rot90(imread([d_dir d_name{1} 'd.jpg']),2);
            
            figure; 
            subplot(2,1,1); imgray(c_im); title(c_name);
            subplot(2,1,2); imgray(d_im); title(d_name{1});
        end
    end
end
%%
create_folder('C:\isbe\nailfold\misc\hilo_mag_ims');
for id = 91:10:171
    
    match = ...
        derm_id_data.people_id == id &...
        derm_id_data.visit == visit &...
        strcmpi(derm_id_data.hand, hand) &...
        derm_id_data.digit == digit;
    
    if any(match)
        
        c_name = derm_id_data.im_names{match};
        d_name = derm_id_data.derm_names{match};
        
        if exist([c_dir c_name '.png'], 'file') && ...
            ~isempty(d_name) && ...
            exist([d_dir d_name{1} 'd.jpg'], 'file')  
        
            c_im = imread([c_dir c_name '.png']);
            d_im = rot90(imread([d_dir d_name{1} 'd.jpg']),2);
            
            imwrite(d_im, ['C:\isbe\nailfold\misc\hilo_mag_ims\' zerostr(id,3) '_lo.png']);
            imwrite(c_im, ['C:\isbe\nailfold\misc\hilo_mag_ims\' zerostr(id,3) '_hi.png']);
            
            [d_roi c_roi] = lo_v_hi_mag_ncm(d_im, c_im);
            write_im_from_colormap([d_roi; c_roi], ['C:\isbe\nailfold\misc\hilo_mag_ims\' zerostr(id,3) '_pair.png'], gray(256));
        end
    end
end

%%
c_im = imread([c_dir '98992c.png']);
d_im = rot90(imread([d_dir '32045d.jpg']),2);

c_roi = c_im(115+(1:480), 175+(1:640));
d_roi = d_im(1335+(1:480), 455+(1:640),:);

figure;
subplot(1,2,1); imgray(c_roi);
subplot(1,2,2); imgray(d_roi);
%%
derm_id_data = u_load('C:\isbe\nailfold\data\rsa_study\data_lists\derm_image_id_data.mat');
load('C:\isbe\nailfold\data\rsa_study\data_lists\miccai_lists.mat', 'miccai_selection');

apex_metrics_dir = 'C:\isbe\nailfold\data\rsa_study\master_set\apex_maps\set12g_half_296655\miccai_maxima\apex_metrics\';

visit2_idx = find(derm_id_data.visit == 2);
num_ims = length(visit2_idx);

visit1_idx = zeros(num_ims,1);
for i_im = 1:num_ims
    
	id = derm_id_data.people_id( visit2_idx(i_im) );
    visit = 1;
    hand = derm_id_data.hand{ visit2_idx(i_im) };
    digit = derm_id_data.digit( visit2_idx(i_im) );
    
    match = find( ...
        derm_id_data.people_id == id &...
        derm_id_data.visit == visit &...
        strcmpi(derm_id_data.hand, hand) &...
        derm_id_data.digit == digit ); 
    
    if ~isempty(match) && ...
            miccai_selection.test(visit2_idx(i_im)) && miccai_selection.test(match) 
            %exist([apex_metrics_dir derm_id_data.im_names{match} '_am.mat'], 'file') && ...
            %exist([apex_metrics_dir derm_id_data.im_names{visit2_idx(i_im)} '_am.mat'], 'file');
                  
        visit1_idx(i_im) = match;
    end
end
missing = ~visit1_idx;
visit1_idx(missing) = [];
visit2_idx(missing) = [];
%%
num_ims = length(visit2_idx);

mean_weighted_width = zeros(num_ims,2);
mean_orientation_entropy = zeros(num_ims,2);
mean_inter_capillary_distance = zeros(num_ims,2);
mean_density = zeros(num_ims,2);

for i_im = 1:num_ims
    visit1_am = u_load([apex_metrics_dir derm_id_data.im_names{visit1_idx(i_im)} '_am.mat']);
    visit2_am = u_load([apex_metrics_dir derm_id_data.im_names{visit2_idx(i_im)} '_am.mat']);
        
    ori_entropy = mb_entropy(visit1_am.distal.orientation_hist,2);
    apex_xy = sortrows(visit1_am.distal.apex_xy);
    inter_d = sqrt(sum(diff(apex_xy).^2,2));
    d2 = sqrt(sum((apex_xy(1,:)-apex_xy(end,:)).^2));
    
    mean_weighted_width(i_im, 1) = naNmedian(visit1_am.distal.mean_weighted_width);
    mean_orientation_entropy(i_im, 1) = naNmedian(ori_entropy);
    mean_inter_capillary_distance(i_im, 1) = median(inter_d);
    mean_density(i_im, 1) = length(visit1_am.distal.mean_weighted_width) / d2;
    
    ori_entropy = mb_entropy(visit2_am.distal.orientation_hist,2);
    apex_xy = sortrows(visit2_am.distal.apex_xy);
    inter_d = sqrt(sum(diff(apex_xy).^2,2));
    d2 = sqrt(sum((apex_xy(1,:)-apex_xy(end,:)).^2));
    
    mean_weighted_width(i_im, 2) = naNmedian(visit2_am.distal.mean_weighted_width);
    mean_orientation_entropy(i_im, 2) = naNmedian(ori_entropy);
    mean_inter_capillary_distance(i_im, 2) = median(inter_d);
    mean_density(i_im, 1) = length(visit2_am.distal.mean_weighted_width) / d2;
end
%
figure;
plot(mean_weighted_width(:,1), mean_weighted_width(:,2), 'x');
axis equal; axis([0 max(mean_weighted_width(:)) 0 max(mean_weighted_width(:))]);
title('Mean weighted width');

figure; plot(mean_orientation_entropy(:,1), mean_orientation_entropy(:,2), 'x');
axis equal; axis([0 max(mean_orientation_entropy(:)) 0 max(mean_orientation_entropy(:))]);
title('Mean orientation entropy');

figure; plot(mean_inter_capillary_distance(:,1), mean_inter_capillary_distance(:,2), 'x');
axis equal; axis([0 max(mean_inter_capillary_distance(:)) 0 max(mean_inter_capillary_distance(:))]);
title('Mean inter-capillary distance');

figure; plot(mean_density(:,1), mean_density(:,2), 'x');
axis equal; axis([0 max(mean_density(:)) 0 max(mean_density(:))]);
title('Mean density');

keep_width = ~any(isnan(mean_weighted_width),2);
rho_width = corr(mean_weighted_width(keep_width,:))

keep_ori = ~any(isnan(mean_orientation_entropy),2);
rho_ori = corr(mean_orientation_entropy(keep_ori,:))

keep_dist = ~any(isnan(mean_inter_capillary_distance),2);
rho_dist = corr(mean_inter_capillary_distance(keep_dist,:))

keep_density = ~any(isnan(mean_density),2);
rho_density = corr(mean_density(keep_density,:))
%%
create_folder('C:\isbe\nailfold\misc\visit_12_ims');
for i_im = 1:num_ims
    im1 = u_load(['C:\isbe\nailfold\data\rsa_study\master_set\images\' derm_id_data.im_names{visit1_idx(i_im)} '.mat']);
    im2 = u_load(['C:\isbe\nailfold\data\rsa_study\master_set\images\' derm_id_data.im_names{visit2_idx(i_im)} '.mat']);
    
    mask1 = u_load(['C:\isbe\nailfold\data\rsa_study\master_set\fov_masks\' derm_id_data.im_names{visit1_idx(i_im)} '_f_mask.mat']);
    mask2 = u_load(['C:\isbe\nailfold\data\rsa_study\master_set\fov_masks\' derm_id_data.im_names{visit2_idx(i_im)} '_f_mask.mat']);
    
    g_lims = [min(min(im1(mask1)), min(im2(mask2))) max(max(im1(mask1)), max(im2(mask2)))];
    
    c1 = size(im1,2);
    c2 = size(im2,2);
    
    if (c1 > c2)
        im2(1:end,c2:c1) = g_lims(2);
    else
        im1(1:end,c1:c2) = g_lims(2);
    end
    
    write_im_from_colormap([im1; im2], ['C:\isbe\nailfold\misc\visit_12_ims\pair' zerostr(i_im,2) '.png'], gray(256), g_lims);
end
%%
c_roi_list = dir('C:\isbe\nailfold\misc\hilo_mag_ims\*_lo.png');

for i_im = 1;%:length(c_roi_list)
    
    c_roi = imread(['C:\isbe\nailfold\misc\hilo_mag_ims\' c_roi_list(i_im).name]);
    c_roi = imresize(c_roi(:,:,2), 0.5, 'lanczos2');
    c_roi_s = imfilter(double(c_roi(201:800,1:800)), fspecial('gaussian', 7, 1), 'replicate');
    figure; imgray(c_roi_s);
    
    detect_capillaries(c_roi_s, 1:7, ...
        'nailfold_mask',            true(size(c_roi_s)),...
        'save_path',                ['C:\isbe\nailfold\misc\hilo_mag_ims\' c_roi_list(i_im).name(1:end-4) '_analysis.mat'],...
        'plot',                     1);
end
%%
detect_capillaries(frame, 1, ...
    'nailfold_mask',            true(size(frame)),...
    'save_path',                'C:\isbe\nailfold\misc\hilo_mag_ims\v_detections',...
    'plot',                     1);



    
    
    


