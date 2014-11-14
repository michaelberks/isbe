% PROBABLY JUNK NOW

clc;
% clear;
close all;
timebar('closeall');

% Make sure the following are created from scratch
clear('gt', 'min_image','mask','du_sum');

imgroot = 'U:\projects\nailfold\capture\2012_10_22\Left.Digit4.x300\';
% imgpath = fullfile(imgroot, 'seq1\preprocessed\registered_g1d\masked');
imgpath = fullfile(imgroot, 'seq2\corrected\registered_g1d\masked');
% imgpath = fullfile(imgroot, '\seq2\corrected\registered_g1d\cropped2');

imgroot = 'U:\projects\nailfold\synthesis';
% d = dir(fullfile(imgroot, '2*T*'));
% imgpath = fullfile(imgroot, d(end).name);
imgpath = fullfile(imgroot, '20130625T101700');

% imgpath = fullfile(imgroot, '');

imgpath = fullfile(imgpath, 'halfsize');

d = dir(fullfile(imgpath,'frame_*.png'));
d = d(2:end);
% d = d(1:10);
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
    min_image = min(min_image, img1);
    min_image = min(min_image, img2);
    
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
%     u = sign(u);  v = sign(v);

    if ~exist('gt','var')
        gt.flowmap = complex(zeros(size(u)), zeros(size(v)));
    end
    
    du = real(gt.flowmap) - u;
    dv = imag(gt.flowmap) - v;
    valid = find(~isnan(du) & ~isnan(dv));

    if ~exist('du_sum','var')
        % Start with empty structures
        n = nan(size(u));
        u_sum = nan(size(u));
        v_sum = nan(size(u)); 
        u_sum_sqr = nan(size(u));
        v_sum_sqr = nan(size(u)); 
        
        % Put in initial values
        n(valid) = 1;
        u_sum(valid) = u(valid);
        v_sum(valid) = v(valid);
        u_sum_sqr(valid) = u(valid).^2;
        v_sum_sqr(valid) = v(valid).^2;
        dn = length(valid);
        du_sum = sum(du(valid));
        du_sum_abs = sum(abs(du(valid)));
        du_sum_sqr = sum(du(valid).^2);
        dv_sum = sum(dv(valid));
        dv_sum_abs = sum(abs(dv(valid)));
        dv_sum_sqr = sum(dv(valid).^2);
        
        u_mat = u;
        v_mat = v;
    else
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
        
        u_mat(:,:,end+1) = u;
        v_mat(:,:,end+1) = v;
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
    rows = 1:size(img1,1);
    cols = 1:size(img1,2);
    
    alpha = 0.3;
%     underlay = 255 * normim(min_image);
%     underlay = 255 * normim(output.img);
    underlay = 255 * double(~mask);
    underlay = repmat(alpha * double(underlay(rows, cols)), [1,1,3]);

    dxy = sqrt(dx.^2 + dy.^2);

    % Overlays should have black background
    bgcolour = [0,0,0];
    dx_rgb = flow_image(dx(rows, cols), bgcolour);
    dy_rgb = flow_image(dy(rows, cols), bgcolour);
    dxy_rgb = flow_image(dxy(rows, cols), bgcolour);

    dx_overlay = (1-alpha) * dx_rgb;
    dy_overlay = (1-alpha) * dy_rgb;
    dxy_overlay = (1-alpha) * dxy_rgb;

    % Redefine colormap so that background is white
    bgcolour = [1,1,1];
    dx_rgb = flow_image(dx(rows, cols), bgcolour);
    dy_rgb = flow_image(dy(rows, cols), bgcolour);
    dxy_rgb = flow_image(dxy(rows, cols), bgcolour);

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
    flowpath = fullfile(flowpath, datestr(now, 'yyyymmddTHHMMSS'));
    
    if ~exist(flowpath, 'dir')
        mkdir(flowpath);
    end
    delete(fullfile(flowpath,'*.png'));
    imwrite(uint8(underlay), fullfile(flowpath,'underlay.png'));
    imwrite(uint8(min_image), fullfile(flowpath,'min_image.png'));
    imwrite(uint8(255*dx_rgb), fullfile(flowpath,'u.png'));
    imwrite(uint8(255*dy_rgb), fullfile(flowpath,'v.png'));
    imwrite(uint8(underlay + dx_overlay), fullfile(flowpath,'u_overlay.png'));
    imwrite(uint8(underlay + dy_overlay), fullfile(flowpath,'v_overlay.png'));
    
    bgcol = 1;
    rgb = complex2rgb(complex(dx, dy), [],[],[], bgcol);
    rgb(isnan(rgb)) = bgcol;
    imwrite(uint8(255*rgb), fullfile(flowpath,'rgb_white.png'));
    bgcol = 0;
    rgb = complex2rgb(complex(dx, dy), [],[],[], bgcol);
    rgb(isnan(rgb)) = bgcol;
    imwrite(uint8(255*rgb), fullfile(flowpath,'rgb_black.png'));
    
    save(fullfile(flowpath, 'uv.mat'), ...
         'u_mat', 'v_mat');
end

if 0
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
