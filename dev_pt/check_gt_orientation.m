function err_stats = check_gt_orientation()

clc;

%% Test accuracy of ground truth orientation of retinograms
% 1. Generate a high resolution test image
% 2. Downsample to a low res image similar to retinogram
% 3. Estimate ground truth orientation using skeletonization
% 4. Compare with 'real' ground truth used to generate the image

imsz = 512;
factor = 8;
nImages = 40;

width_vec = 2.^linspace(2, 4, 16);
nWidths = length(width_vec);

error_vec = [];
plot_vec = zeros(nWidths,4);

for iWidth = 1:nWidths
    width = width_vec(iWidth);
    
    for iImage = 1:nImages
%% 1. Generate a high resolution test image
        % 1. Ask user for (x,y) inputs in figure

        % nknots = 7;
        % [x,y] = ginput(nknots);

        % 2. Generate a smooth b-spline using given points as knot points
        imghi = zeros(imsz,imsz);
        % y1 = ceil(imsz/2 + rand*imsz/4);
        % y2 = ceil(imsz/2 + rand*imsz/4);
        % for x = 1:imsz
        %     y = round(y1 + (y2-y1)/imsz * x);
        %     imghi(y,x) = 1;
        % end
        % ori_gt = -atan2(y2-y1, imsz) * 180/pi;
        % disp(ori_gt);
        [imghi, lbl, lblc, ori_gt] = create_sin_curve(width, 255, ...
                                                      imsz/4 + rand*imsz/2, ...
                                                      rand*180/pi, 1, ...
                                                      imsz, imsz, ...
                                                      imsz/2, imsz/2);

%% 2. Downsample to a low res image similar to retinogram
        imglo = imghi(1:factor:end, 1:factor:end);
        ori_gt = ori_gt(1:factor:end, 1:factor:end);
        imshow(imglo*255);

%% 3. Estimate ground truth orientation using skeletonization
        % 1. Skeletonize the line image
        % 2. Create locally linear fits along the line
        ori_gt = complex(cosd(2*ori_gt), sind(2*ori_gt));
        ori_gt = ori_gt(imglo ~= 0);
        ori_est = get_gt_ori(imglo);
        ori_est = ori_est(imglo ~= 0);

%% 4. Compare with 'real' ground truth used to generate the image
        ori_err = ori_error(ori_gt, ori_est);
        error_vec = [error_vec; ori_err(~isnan(ori_err))];
    end

    err_stats = ori_error_stats(error_vec);
    plot_vec(iWidth,:) = [width err_stats.abs_mean ...
                                err_stats.abs_range];
end

figure(1); 
    imshow(imglo);
    
figure(2); clf; hold on;
    plot(plot_vec(:,1), plot_vec(:,2),'-');
    plot(plot_vec(:,1), plot_vec(:,3),':');
%     plot(plot_vec(:,1), plot_vec(:,4),':');


function gt_ori = get_gt_ori(gt)

% Model parameters
win_size = 11;

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
        uni_y = sparse(xi, 1, yi, win_size, 1) ./ ...
                sparse(xi, 1, 1, win_size, 1);
        uni_y = full(uni_y(uni_x));
    else
        uni_x = sparse(yi, 1, xi, win_size, 1) ./ ...
                sparse(yi, 1, 1, win_size, 1);
        uni_x = full(uni_x(uni_y));
    end

    uu = mean(diff(uni_x));
    vv = -mean(diff(uni_y));
    dd = sqrt(uu^2 + vv^2);
    gts_ori(c_y(ii), c_x(ii)) = complex(uu / dd, vv / dd);
end

% Assign orientation of nearest centreline to all line points
a_u = griddata(c_x, c_y, real(gts_ori(gts)),a_x, a_y, 'nearest');
a_v = griddata(c_x, c_y, imag(gts_ori(gts)),a_x, a_y, 'nearest');

% Note double angle in complex form
gt_ori = zeros(size(gt));
gt_ori((gt~=0)) = (complex(a_u, a_v).^2) ./ (a_u.^2 + a_v.^2);

