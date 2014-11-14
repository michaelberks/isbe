win_size = 11;
mkdir C:\isbe\asymmetry_project\data\retinograms\DRIVE\training\orientations
mkdir C:\isbe\asymmetry_project\data\retinograms\DRIVE\test\orientations

for jj = 1:40     
    if jj < 21
        data = 'test';
    else
        data = 'training';
    end
    
    gt = imread(['C:\isbe\asymmetry_project\data\retinograms\DRIVE\' data '\1st_manual\' zerostr(jj,2) '_manual1.gif']);
    gt = gt > 0;
    gts = bwmorph(gt, 'skel', 'inf');


%     figure; 
%     subplot(2,2,1); imagesc(gt); axis image;
%     subplot(2,2,2); imagesc(gts); axis image;

    [c_y c_x] = find(gts);
    [a_y a_x] = find(gt);

    gts_ori = zeros(size(gts));
    
    %
    for ii = 1:size(c_x,1)

        local_win = sample_window(gts, win_size, c_y(ii), c_x(ii), 0);

        [yi xi] = find(bwselect(local_win, (win_size+1)/2, (win_size+1)/2, 8));
        uni_x = unique(xi);
        uni_y = unique(yi);

        if length(uni_x) > length(uni_y)
            uni_y = sparse(xi, 1, yi, win_size, 1) ./ sparse(xi, 1, 1, win_size, 1);
            uni_y = full(uni_y(uni_x));
        else
            uni_x = sparse(yi, 1, xi, win_size, 1) ./ sparse(yi, 1, 1, win_size, 1);
            uni_x = full(uni_x(uni_y));
        end

        uu = mean(diff(uni_x));
        vv = -mean(diff(uni_y));
        dd = sqrt(uu^2 + vv^2);
        gts_ori(c_y(ii), c_x(ii)) = complex(uu / dd, vv / dd);
    end
    
    a_u = griddata(c_x, c_y, real(gts_ori(gts)),a_x, a_y, 'nearest');
    a_v = griddata(c_x, c_y, imag(gts_ori(gts)),a_x, a_y, 'nearest');
    
    dd = a_u.^2 + a_v.^2;

    gt_ori = zeros(size(gt));
    gt_ori(gt) = (complex(a_u, a_v).^2) ./ dd;
    
    save(['C:\isbe\asymmetry_project\data\retinograms\DRIVE\' data '\orientations\' zerostr(jj,2) '_ori1.mat'],...
        'gt_ori');
    %figure; image(complex2rgb(gt_ori.^2)); axis image;
end
%%
spacing = 4;
compute_k_maps_batch('2004_screening_processed/abnormals', 'rf', 'spacing', spacing, 'num_jobs', 275, 'task_id', 53);
load C:\isbe\asymmetry_project\data\k_stellate_maps\rf\2004_screening_processed\abnormals\024RCC_mask
f1 = load_uint8('C:\isbe\asymmetry_project\data\k_stellate_maps\rf\2004_screening_processed\abnormals\024RCC_f1');
f2 = load_uint8('C:\isbe\asymmetry_project\data\k_stellate_maps\rf\2004_screening_processed\abnormals\024RCC_f2');
for lev = 1:5
    f1_map = zeros(size(mask));
    f1_map(mask) = f1(:,lev);
    
    f2_map = zeros(size(mask));
    f2_map(mask) = f2(:,lev);
    
    figure; 
    subplot(1,2,1); imagesc(f1_map(1:spacing:end,1:spacing:end)); axis image;
    subplot(1,2,2); imagesc(f2_map(1:spacing:end,1:spacing:end)); axis image;
end
%%
map_dir = 'C:\isbe\asymmetry_project\data\orientation_maps\rf\2004_screening_processed\abnormals\';
mask_dir = 'C:\isbe\asymmetry_project\data\masks\2004_screening\abnormals\';
pec_dir = 'C:\isbe\asymmetry_project\data\masks_pectoral\2004_screening_processed\abnormals\';

map_list = dir([map_dir '*ML*.mat']);
mask_names = match_mammo_names(mask_dir, get_mammo_info(map_list));

for ii = 1:20
    ori_map = angle(load_uint8([map_dir map_list(ii).name]));
    mask = u_load([mask_dir mask_names{ii}]);
    pec = u_load([pec_dir mask_names{ii}]);
    
    [discard_mask] = discard_orientations(ori_map, 'mask', mask & ~pec, 'edge_size', 100);
    figure; imagesc(discard_mask); axis image;
    
    clear discard_mask mask ori_map;
end
%%
mam_dir = 'C:\isbe\asymmetry_project\data\mammograms\2004_screening\normals\';
mam_list = dir([mam_dir '*ML*.mat']);
mask_dir = 'C:\isbe\asymmetry_project\data\masks\2004_screening_processed\normals\';
pec_dir = 'C:\isbe\asymmetry_project\data\masks_pectoral\2004_screening_processed\normals\';
mask_names = match_mammo_names(mask_dir, get_mammo_info(mam_list));
%%
for ii = 36:length(mam_list)
    mam = u_load([mam_dir mam_list(ii).name]);
    mam_small = imresize(mam, [1024 nan]);
    
    pectoral_mask = pectoral_user(mam_small);
    pectoral_mask = imresize(pectoral_mask, size(mam));
    
    save([pec_dir mask_names{ii}], 'pectoral_mask');
    
    answer = questdlg(...
            'Continue to next image','Mammogram Complete','Yes', 'No', 'Yes');
    
    if strcmpi(answer, 'no')
        display(ii);
        break;
    end
end
%%
% mam_dir = 'C:\isbe\asymmetry_project\data\mammograms\2004_screening\abnormals\';
% mam_names = get_mammo_info(dir([mam_dir '*.mat']));
mam_names = {'024RCC', '024RML'};
%map_dir = 'Z:\data\k_stellate_maps\rf\2004_screening_processed\abnormals\';
map_dir = 'C:\isbe\asymmetry_project\data\k_stellate_maps\rf\2004_screening_processed\abnormals\';
spacing = 2;
for ii = 1:2
    
    mask = u_load([map_dir mam_names{ii} '_mask.mat']);
    f1 = load_uint8([map_dir mam_names{ii} '_f1.mat']);
    f2 = load_uint8([map_dir mam_names{ii} '_f2.mat']);
    
    for lev = 1:5
        f1_map = zeros(size(mask));
        f2_map = zeros(size(mask));
        
        f1_map(mask) = f1(:,lev);
        f2_map(mask) = f2(:,lev);
        
        figure;
        subplot(1,2,1); imagesc(f1_map(1:spacing:end, 1:spacing:end)); axis image;
        subplot(1,2,2); imagesc(f2_map(1:spacing:end, 1:spacing:end)); axis image;
    end
end
%%
num_pts = zeros(20,1);
num_pt = zeros(20,1);
for jj = 21:40     
    if jj < 21
        data = 'test';
    else
        data = 'training';
    end
    
    gt = imread(['C:\isbe\asymmetry_project\data\retinograms\DRIVE\' data '\1st_manual\' zerostr(jj,2) '_manual1.gif']);
    gt = gt > 0;
    gts = bwmorph(gt, 'skel', 'inf');

    num_pts(jj-20) = sum(gts(:));
    num_pt(jj-20) = sum(gt(:));
end
        
        
    
    