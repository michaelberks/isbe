% PROBABLY JUNK NOW

clear;

% datroot = 'U:\projects\nailfold\synthesis\20130320T162038\';
% gt = load(fullfile(datroot, '_ground_truth.mat'));
% 
% u_gt = real(gt.flowmap);
% v_gt = imag(gt.flowmap);
% 
% datpath = fullfile(datroot, 'flow');
% d = dir(fullfile(datpath, '*_stats.mat'));
% datfile = fullfile(datpath, d(end).name);
% 
% load(datfile);
% 
% cmap = [ zeros(1,3); redgreen(254); zeros(1,3) ];
% 
% d = dir(fullfile(datroot, 'frame_*.png'));
% img1 = imread(fullfile(datroot, d(241).name));
% 
% for i = 1:size(stats.v_mat,3)
%     img2 = imread(fullfile(datroot, d(i+241).name));
%     v = stats.v_mat(:,:,i);
%     v_img = 1 + 253*(0.5 + v/4);
%     
%     c_inds = ~isnan(min(v,[],1));
%     r_inds = ~isnan(min(v,[],2));
% 
%     figure(1); clf; colormap(cmap);
%         image(uint8(v_img(r_inds, c_inds))); axis('image');
%         
%     figure(2); clf; colormap(gray(256));
%         subplot(1,2,1); imagesc(img1(r_inds, c_inds)); axis('image');
%         subplot(1,2,2); imagesc(img2(r_inds, c_inds)); axis('image');
%     
%     v_sub = v(r_inds, c_inds);
%     v_sub_gt = v_gt(r_inds, c_inds);
%         
%     pause;
%     
%     img1 = img2;
% end

clc;
% clear;
close all;
timebar('closeall');

% imgroot = 'U:\projects\nailfold\capture\2012_10_22\Left\Digit4\x300\';
% imgpath = fullfile(imgroot, 'seq2\corrected\registered_g1d');
% 
% mosaic = double(imread(fullfile(imgpath,'mosaic.png')));
% mosaic(mosaic==0) = nan;
% 
% mosaic_n = 255*normim(mosaic);
% mosaic_n(isnan(mosaic_n)) = 0;
% underlay = repmat(double(mosaic_n), [1,1,3]);
% overlay = zeros(size(mosaic));
% 
% n_sigma = 3;
% line_threshold = 3;
% for i = n_sigma%:-1:1
%     sigma = 2^i;
% %     [line_strength, line_orientation] = ...
% %         gaussian_1st_derivative_gradient(mosaic, sigma);
% %     line_strength = abs(line_strength);
% 
%     [line_strength, line_orientation] = ...
%         gaussian_2nd_derivative_line(mosaic, sigma);
% 
%     line_nms = mb_non_maximal_supp(line_strength, line_orientation);
%     line_nms(line_nms < line_threshold) = 0;
%     
%     overlay(abs(line_nms) > 0) = i;
% end
% overlay = 255 * ind2rgb(overlay+1, [0,0,0; jet(n_sigma)]);
% 
% alpha = 0.8;
% img = alpha * underlay + (1-alpha) * overlay;
% figure('windowstyle','docked'); clf;
%     imagesc(uint8(img)); axis('image');
% figure('windowstyle','docked'); clf;
%     imagesc(uint8(alpha*underlay)); axis('image');
%     
% return
% 
% hw = 15; % half width
% fw = 2*hw + 1; % full width
% 
% img0 = zeros(fw,fw);
% % img0(:, hw:end) = 1; % X only
% img0(hw:end, :) = 1; % Y only
% % img0(hw:end, hw:end) = 1; % X and Y
% 
% % f = isbe_fspecial('gaussian', [9,9], 3);
% f = ones(9,9); f = f / sum(abs(f(:)));
% img0 = conv2(img0, f, 'valid');
% 
% img1 = img0(2:end, 2:end);
% img2 = img0(1:end-1, 1:end-1);
% 
% figure(1); clf; hold on; colormap(gray(256));
%     imagesc(img1); axis('image','ij');
% figure(2); clf; hold on; colormap(gray(256));
%     imagesc(img2); axis('image','ij');
% 
% Ix = conv2(img1, 0.5*[0 0 0; 1 0 -1; 0 0 0], 'valid');
% Iy = conv2(img1, 0.5*[0 0 0; 1 0 -1; 0 0 0]', 'valid');
% It = (img2 - img1); It = It(2:end-1, 2:end-1);
% 
% npts = numel(Ix);
% z = zeros(npts,1);
% 
% figure(3); clf; hold on;
%     plot3([z Ix(:)]', [z Iy(:)]', [z It(:)]', 'b-');
%     xyzlabel('Ix','Iy','It');
%     axis('equal');
% 
% return

% Make sure the following are created from scratch
clear('min_image','mask','du_sum');

imgroot = 'U:\projects\nailfold\capture\2012_10_22\Left.Digit4.x300\';
% imgpath = fullfile(imgroot, 'seq1\preprocessed\registered_g1d\masked');
imgpath = fullfile(imgroot, 'seq2\corrected\registered_g1d\masked');

% imgroot = 'U:\tmp\ncm\synthesis';
% imgpath = fullfile(imgroot, '');

imgpath = fullfile(imgpath, 'halfsize');

d = dir(fullfile(imgpath,'frame_*.png'));
d = d(1:20);
% d = d(1:3:end);

if strcmp(username(), 'ptresadern')
    tb = timebar('title', 'Computing flow', ...
                 'limit', length(d)-1);
end

if 1
    gt_filename = fullfile(imgpath,'ground_truth.mat');
    if exist(gt_filename, 'file')
        gt = load(gt_filename);
    end
end
    
for i = 2:length(d)
    img1 = mean(imread(fullfile(imgpath, d(i-1).name)), 3);
    img2 = mean(imread(fullfile(imgpath, d(i).name)), 3);
    
    if ~exist('min_image','var')
        min_image = inf(size(img1));
    end
    min_image = min(min(min_image, img1), img2);
    
    patch_hw = 2;
    output = pt_flow_derivs(img1, img2, [], [], ...
                            patch_hw);
   
    if ~exist('mask','var')
        mask = false(size(img1));
    end
    mask = mask | ...
           ( (output.u ~= 0) & ~isnan(output.u) | ...
             (output.v ~= 0) & ~isnan(output.v) );
    
    u = output.u; v = output.v;
    u = sign(u);  v = sign(v);

    du = real(gt.flowmap) - u;
    dv = imag(gt.flowmap) - v;
    valid = find(~isnan(du) & ~isnan(dv));

    if ~exist('du_sum','var')
        n = nan(size(u));
        n(valid) = 1;
        
        u_sum = nan(size(u));
        u_sum(valid) = u(valid);
        v_sum = nan(size(u)); 
        v_sum(valid) = v(valid);
        
        u_sum_sqr = nan(size(u));
        u_sum_sqr(valid) = u(valid).^2;
        v_sum_sqr = nan(size(u)); 
        v_sum_sqr(valid) = v(valid).^2;

        dn = length(valid);
        du_sum = sum(du(valid));
        du_sum_abs = sum(abs(du(valid)));
        du_sum_sqr = sum(du(valid).^2);
        dv_sum = sum(dv(valid));
        dv_sum_abs = sum(abs(dv(valid)));
        dv_sum_sqr = sum(dv(valid).^2);
    else
        if 1
            % Add 'em up
            n(valid) = n(valid) + 1;

            u_sum(valid) = u_sum(valid) + u(valid);
            v_sum(valid) = v_sum(valid) + v(valid);
            
            u_sum_sqr(valid) = u_sum_sqr(valid) + u(valid).^2;
            v_sum_sqr(valid) = v_sum_sqr(valid) + v(valid).^2;

            dn = dn + length(valid);
            du_sum = du_sum + sum(du(valid));
            du_sum_abs = du_sum_abs + sum(abs(du(valid)));
            du_sum_sqr = du_sum_sqr + sum(du(valid).^2);
            dv_sum = dv_sum + sum(dv(valid));
            dv_sum_abs = dv_sum_abs + sum(abs(dv(valid)));
            dv_sum_sqr = dv_sum_sqr + sum(dv(valid).^2);
        else
            % Take the value with the higher absolute value
            inds = find(abs(u(valid)) > abs(dx(valid)));
            dx(valid(inds)) = u(valid(inds));
            inds = find(abs(v(valid)) > abs(dy(valid)));
            dy(valid(inds)) = v(valid(inds));
        end
    end
    
    if exist('tb', 'var')
        timebar(tb, 'advance');
    end
end
if exist('tb', 'var')
    timebar(tb, 'close');
end

dx = u_sum ./ n;
dy = v_sum ./ n;

dx_var = (u_sum_sqr./n) - dx.^2;
dy_var = (v_sum_sqr./n) - dy.^2;

du_mean = du_sum / dn;
du_mean_abs = du_sum_abs / dn;
du_mean_sqr = du_sum_sqr / dn;
du_rms = sqrt(du_mean_sqr);
du_var = du_mean_sqr - du_mean^2;
du_std = sqrt(du_var);

dv_mean = dv_sum / dn;
dv_mean_abs = dv_sum_abs / dn;
dv_mean_sqr = dv_sum_sqr / dn;
dv_rms = sqrt(dv_mean_sqr);
dv_var = dv_mean_sqr - dv_mean^2;
dv_std = sqrt(dv_var);


%% Show independent distributions of estimated motions
if 0
    interior0 = output.interior;
    [m,n] = size(interior0);
    interior = false([m,n]);
    for y = 2:m-1
        for x = 2:n-1
            interior(y,x) = all(all( interior0((y-1):(y+1), (x-1):(x+1)) ));
        end
    end

    figure; maximize();
        imagesc(interior); axis('image');
    figure(); 
        hist(dx(interior), linspace(-2,2,51));
    figure(); 
        hist(dy(interior), linspace(-2,2,51));

    return
end

%% Show joint distribution of estimated motions 
if 0
    n_disp = 101;
    max_disp = 0.5;
    dx_rng = linspace(-max_disp, max_disp, n_disp);
    dy_rng = linspace(-max_disp, max_disp, n_disp);
    hough = hist3([output.u(:) output.v(:)], 'edges', {dx_rng, dy_rng});

    figure; clf; maximize();
        imagesc(dx_rng, dy_rng, hough);
end

%% Show quiver map
if 0
    figure(1); clf; colormap(gray(256)); maximize();
        subplot(2,2,1); imagesc(img1); axis('image');
        subplot(2,2,2); imagesc(img2); axis('image');
        subplot(2,2,3); hold on;
        imagesc(img1);
        quiver(output.x_vec, output.y_vec, output.u_vec, output.v_vec, 5);
        axis('image','ij');
    linkaxes([subplot(2,2,1), subplot(2,2,2), subplot(2,2,3)]);
end

%% Show gradient images
if 0
    figure(2); clf; colormap(gray(256)); maximize();
        subplot(2,2,1); imagesc(output.Ix); axis('image');
        subplot(2,2,2); imagesc(output.Iy); axis('image');
        subplot(2,2,3); imagesc(output.It); axis('image');
    linkaxes([subplot(2,2,1), subplot(2,2,2), subplot(2,2,3)]);
end

%% Colormapped flow images
if 1
%     rows = ~all(isnan(dx),2); 
%     cols = ~all(isnan(dx),1);
    rows = 1:size(img1,1);
    cols = 1:size(img1,2);
    
    alpha = 0.3;
%     underlay = 255 * normim(min_image);
%     underlay = 255 * normim(output.img);
    underlay = 255 * double(~mask);
    underlay = repmat(alpha * double(underlay(rows, cols)), [1,1,3]);

    cmap = [0 0 0; redgreen(255)];

    dx2 = normim(dx(rows, cols), 'stretch_fixed');
    dx_rgb = 255 * ind2rgb(uint8(ceil(1+254*dx2)), cmap);
    dx_overlay = (1-alpha) * dx_rgb;

    dy2 = normim(dy(rows, cols), 'stretch_fixed');
    dy_rgb = 255 * ind2rgb(uint8(ceil(1+254*dy2)), cmap);
    dy_overlay = (1-alpha) * dy_rgb;

    dxy = sqrt(dx.^2 + dy.^2);
    dxy2 = normim(dxy(rows, cols), 'stretch_fixed');
    dxy_rgb = 255 * ind2rgb(uint8(ceil(1+254*dxy2)), cmap);
    dxy_overlay = (1-alpha) * dxy_rgb;

    % Redefine colormap so that background is white
    cmap = [1 1 1; redgreen(255)];
    dx_rgb = 255 * ind2rgb(uint8(ceil(1+254*dx2)), cmap);
    dy_rgb = 255 * ind2rgb(uint8(ceil(1+254*dy2)), cmap);
    dxy_rgb = 255 * ind2rgb(uint8(ceil(1+254*dxy2)), cmap);

%     figure; clf; colormap(gray(256)); maximize();
%         image(uint8(underlay + dx_overlay)); axis('image');
%     figure; clf; maximize();
%         output.var_u(output.var_u > 1.0) = 0;
%         imagesc(log10(output.var_u)); axis('image');

%     figure; clf; colormap(gray(256)); maximize();
%         image(uint8(underlay + dy_overlay)); axis('image');
%     figure; clf; maximize();
%         output.var_v(output.var_v > 1.0) = 0;
%         imagesc(log10(output.var_v)); axis('image');
        
%     figure; clf; colormap(gray(256)); maximize();
%         image(uint8(underlay + dxy_overlay)); axis('image');
        
%     flowpath = 'M:\nailfold\weekly_presentations\2012-12-12\flow';
    flowpath = fullfile(imgpath,'flow');
    
    if ~exist(flowpath, 'dir')
        mkdir(flowpath);
    end
    delete(fullfile(flowpath,'*.png'));
    imwrite(uint8(underlay), fullfile(flowpath,'underlay.png'));
    imwrite(uint8(min_image), fullfile(flowpath,'min_image.png'));
    imwrite(uint8(dx_rgb), fullfile(flowpath,'u.png'));
    imwrite(uint8(dy_rgb), fullfile(flowpath,'v.png'));
    imwrite(uint8(underlay + dx_overlay), fullfile(flowpath,'u_overlay.png'));
    imwrite(uint8(underlay + dy_overlay), fullfile(flowpath,'v_overlay.png'));
end

if 1
    gt_filename = fullfile(imgpath,'ground_truth.mat');
    if exist(gt_filename, 'file')
        gt = load(gt_filename);
    end
    
    img = [real(gt.flowmap) dx imag(gt.flowmap) dy];
    img = 255 * normim(img, 'stretch_fixed');
    figure; colormap([0 0 0; redgreen(255)]);
        image(uint8(img)); axis('image');
        
    err_x = dx - real(gt.flowmap);
    mae_u = mean(abs(err_x(~isnan(err_x))));
    err_y = dy - imag(gt.flowmap);
    mae_v = mean(abs(err_y(~isnan(err_y))));
    
    mean_u_var = mean(dx_var(~isnan(dx_var)));
    mean_v_var = mean(dy_var(~isnan(dy_var)));
    
    mean_u_std = mean(sqrt(dx_var(~isnan(dx_var))));
    mean_v_std = mean(sqrt(dy_var(~isnan(dy_var))));
    
    disp([mae_u mean_u_var mean_u_std; 
          mae_v mean_v_var mean_v_std]);
        
    datestring = datestr(now,'yyyymmddTHHMMSS');

    filename = sprintf('%s_stats.mat', datestring);
    parameters = gt.parameters;
    
    if exist('flowsubdir','var')
        flowpath = fullfile(flowpath, flowsubdir);
    end
    
    save(fullfile(flowpath,filename), ...
         'parameters', 'mae_u', 'mae_v', 'dx_var', 'dy_var');
    
    filename = sprintf('%s_stats.txt', datestring);
    fid = fopen(fullfile(flowpath,filename), 'w');
    if (fid ~= 0)
     	fprintf(fid, '%s\n', evalc('gt.parameters'));
        fprintf(fid, 'MAE(u) = %5.3f\n', mae_u);
        fprintf(fid, 'MAE(v) = %5.3f\n', mae_v);
        fprintf(fid, 'mean(var(u)) = %5.3f\n', mean_u_var);
        fprintf(fid, 'mean(var(v)) = %5.3f\n', mean_v_var);
        fprintf(fid, 'mean(std(u)) = %5.3f\n', mean_u_std);
        fprintf(fid, 'mean(std(v)) = %5.3f\n', mean_v_std);
        fprintf(fid, '\n');
        fprintf(fid, 'du_mean = %5.3f\n', du_mean);
        fprintf(fid, 'du_mean_abs = %5.3f\n', du_mean_abs);
        fprintf(fid, 'du_mean_sqr = %5.3f\n', du_mean_sqr);
        fprintf(fid, 'du_rms = %5.3f\n', du_rms);
        fprintf(fid, 'du_var = %5.3f\n', du_var);
        fprintf(fid, 'du_std = %5.3f\n', du_std);
        fprintf(fid, '\n');
        fprintf(fid, 'dv_mean = %5.3f\n', dv_mean);
        fprintf(fid, 'dv_mean_abs = %5.3f\n', dv_mean_abs);
        fprintf(fid, 'dv_mean_sqr = %5.3f\n', dv_mean_sqr);
        fprintf(fid, 'dv_rms = %5.3f\n', dv_rms);
        fprintf(fid, 'dv_var = %5.3f\n', dv_var);
        fprintf(fid, 'dv_std = %5.3f\n', dv_std);
        fclose(fid);
    end
    return;
    
    flowmap = gt.flowmap;
    u_gt = 255 * normim(real(flowmap), 'stretch_fixed');
    v_gt = 255 * normim(imag(flowmap), 'stretch_fixed');
    
    figure; colormap(redgreen(256));
        sp = [1,4];
        subplot(sp(1),sp(2),1); image(uint8(u_gt)); axis('image');
        subplot(sp(1),sp(2),2); image(uint8(dx_overlay)); axis('image');
        subplot(sp(1),sp(2),3); image(uint8(v_gt)); axis('image');
        subplot(sp(1),sp(2),4); image(uint8(dy_overlay)); axis('image');
end
