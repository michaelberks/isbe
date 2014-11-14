%compute GT orientations for the DRIVE data
win_size = 11;
win_size_m = 7; %use smaller win size for mixed orientations
mkdir C:\isbe\asymmetry_project\data\retinograms\STARE\training\orientations
mkdir C:\isbe\asymmetry_project\data\retinograms\STARE\training\mixed_orientations

for jj = 1:20
    
    %Load in ground truth and skeletonise
    gt = u_load(['C:\isbe\asymmetry_project\data\retinograms\STARE\training\vessel_masks\' zerostr(jj,2) '_training_v_mask.mat']);
    gts = bwmorph(gt, 'skel', 'inf');

    %Extract x,y coords of vessels + vessel centres
    [c_y c_x] = find(gts);
    [a_y a_x] = find(gt);

    %Create storage for the ground truth orientations
    gts_ori = zeros(size(gts));
    
    num_pts = size(c_x,1);
    
    %Loop through each skeleton point
    for ii = 1:num_pts

        %Sample local window from skeleton map
        local_win = sample_window(gts, win_size, c_y(ii), c_x(ii), 0);

        %Get all points connected to the centre
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
    
    %Now loop through the points again to create mixed orientations
    mixed_oris = cell(num_pts,1);
    for ii = 1:num_pts

        %Sample local window from skeleton map
        local_win = sample_window(gts, win_size_m, c_y(ii), c_x(ii), 0);
        local_ori = sample_window(gt_ori, win_size_m, c_y(ii), c_x(ii), 0);
        local_connections = bwselect(local_win, (win_size_m+1)/2, (win_size_m+1)/2, 8);

        mixed_oris{ii} = local_ori(local_connections);
    end
    a_idx = griddata(c_x, c_y, 1:num_pts,a_x, a_y, 'nearest');
    
    mixed_idx = zeros(size(gts));
    mixed_idx(gt) = a_idx;
    
    save(['C:\isbe\asymmetry_project\data\retinograms\STARE\training\mixed_orientations\' zerostr(jj,2) '_training_mixed_ori.mat'],...
        'mixed_oris', 'mixed_idx');
    save(['C:\isbe\asymmetry_project\data\retinograms\STARE\training\orientations\' zerostr(jj,2) '_training_ori.mat'],...
        'gt_ori');
    figure; image(complex2rgb(gt_ori)); axis image;
end
%%
mkdir C:\isbe\asymmetry_project\data\retinograms\STARE\training\width_maps

for jj = 1:20
    
    %Load in ground truth and skeletonise
    gt = u_load(['C:\isbe\asymmetry_project\data\retinograms\STARE\training\vessel_masks\' zerostr(jj,2) '_training_v_mask.mat']);
    gts = bwmorph(gt, 'skel', 'inf');

    %Extract x,y coords of vessels + vessel centres
    [c_y c_x] = find(gts);
    [a_y a_x] = find(gt);

    %Create storage for the ground truth orientations
    gts_width = double(bwdist(~gt));
    
    a_w = griddata(c_x, c_y, gts_width(gts), a_x, a_y, 'nearest');

    gt_width = zeros(size(gt));
    gt_width(gt) = a_w;

    save(['C:\isbe\asymmetry_project\data\retinograms\STARE\training\width_maps\' zerostr(jj,2) '_training_width.mat'],...
        'gt_width');
    figure; image(gt_width); axis image;
end