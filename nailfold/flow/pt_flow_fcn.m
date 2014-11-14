% PROBABLY JUNK NOW

function pt_flow_fcn(imgpath, flowsubdir)

if (nargin==0 && nargout==0), test_script(); return; end

d = dir(fullfile(imgpath,'frame_*.png'));
d = d(240:739);
% d = d(1:100);

gt_filename = fullfile(imgpath,'_ground_truth.mat');
if exist(gt_filename, 'file')
    gt = load(gt_filename);
end

if strcmp(username(), 'ptresadern')
    tb = timebar('title', 'Computing flow', ...
                 'limit', length(d)-1);
end

img1 = mean(imread(fullfile(imgpath, d(1).name)), 3);
min_image = img1;
mask = false(size(img1));

Ix = []; 
Iy = [];

% Initialize counters
z = nan(size(img1));

nu = z;
u_sum = z;
u_sum_sqr = z;
du_sum = z;
du_sum_abs = z;
du_sum_sqr = z;

nv = z;
v_sum = z;
v_sum_sqr = z;
dv_sum = z;
dv_sum_abs = z;
dv_sum_sqr = z;

u_mat = zeros([size(img1) length(d)-1]);
v_mat = zeros([size(img1) length(d)-1]);

for i = 2:length(d)
    img2 = mean(imread(fullfile(imgpath, d(i).name)), 3);
    min_image = min(min_image, img2);
    
    patch_hw = 3;
    if 0
        tic;
        output = pt_flow_derivs(img1, img2, Ix, Iy, ...
                                patch_hw);
        toc;
    else
        output = pt_flow_derivs_joint(img1, img2, Ix, Iy, ...
                                      1, ... % patch halfwidth
                                      1, ... % step
                                      0.01, ... % lambda
                                      0.0 ... % tau
                                     );
    end
    
    mask = mask | ( (output.u ~= 0) & ~isnan(output.u) | ...
                    (output.v ~= 0) & ~isnan(output.v) );
    
    u = output.u;
    v = output.v;
    
%     figure(1);
%     subplot(1,2,2); image(flow_image(u)); axis('image');
%     return;
    
    du = real(gt.flowmap) - u;
    dv = imag(gt.flowmap) - v;

    u_valid = (~isnan(du) & (gt.mask>0));
    v_valid = (~isnan(dv) & (gt.mask>0));

    % Add 'em up
    inds = (u_valid & ~isnan(u_sum));
        nu(inds) = nu(inds) + 1;
        u_sum(inds) = u_sum(inds) + u(inds);
        u_sum_sqr(inds) = u_sum_sqr(inds) + u(inds).^2;
        du_sum(inds) = du_sum(inds) + du(inds);
        du_sum_abs(inds) = du_sum_abs(inds) + abs(du(inds));
        du_sum_sqr(inds) = du_sum_sqr(inds) +     du(inds).^2;
    inds = (u_valid & isnan(u_sum));
        nu(inds) = 1;
        u_sum(inds) = u(inds);
        u_sum_sqr(inds) = u(inds).^2;
        du_sum(inds) = du(inds);
        du_sum_abs(inds) = abs(du(inds));
        du_sum_sqr(inds) =     du(inds).^2;

    inds = (v_valid & ~isnan(v_sum));
        nv(inds) = nv(inds) + 1;
        v_sum(inds) = v_sum(inds) + v(inds);
        v_sum_sqr(inds) = v_sum_sqr(inds) + v(inds).^2;
        dv_sum(inds) = dv_sum(inds) + dv(inds);
        dv_sum_abs(inds) = dv_sum_abs(inds) + abs(dv(inds));
        dv_sum_sqr(inds) = dv_sum_sqr(inds) +     dv(inds).^2;
    inds = (v_valid & isnan(v_sum));
        nv(inds) = 1;
        v_sum(inds) = v(inds);
        v_sum_sqr(inds) = v(inds).^2;
        dv_sum(inds) = dv(inds);
        dv_sum_abs(inds) = abs(dv(inds));
        dv_sum_sqr(inds) =     dv(inds).^2;

    u_mat(:,:,i-1) = u;
    v_mat(:,:,i-1) = v;

    img1 = img2;
    Ix = output.Ix;
    Iy = output.Iy;
    
    if exist('tb', 'var')
        timebar(tb, 'advance');
    end
end
if exist('tb', 'var')
    timebar(tb, 'close');
end

stats.u_mat = u_mat;
stats.u_mean = u_sum ./ nu;
stats.u_var = (u_sum_sqr ./ nu) - (stats.u_mean).^2;
stats.du_mean = du_sum ./ nu;
stats.du_mean_abs = du_sum_abs ./ nu;
stats.du_mean_sqr = du_sum_sqr ./ nu;
stats.du_rms = sqrt(stats.du_mean_sqr);
stats.du_var = stats.du_mean_sqr - (stats.du_mean).^2;
stats.du_std = sqrt(stats.du_var);

stats.v_mat = v_mat;
stats.v_mean = v_sum ./ nv;
stats.v_var = (v_sum_sqr ./ nv) - (stats.v_mean).^2;
stats.dv_mean = dv_sum ./ nv;
stats.dv_mean_abs = dv_sum_abs ./ nv;
stats.dv_mean_sqr = dv_sum_sqr ./ nv;
stats.dv_rms = sqrt(stats.dv_mean_sqr);
stats.dv_var = stats.dv_mean_sqr - (stats.dv_mean).^2;
stats.dv_std = sqrt(stats.dv_var);

flowpath = fullfile(imgpath,'flow');
if ~exist(flowpath, 'dir')
    mkdir(flowpath);
end

% show_interior(output, dx, dy);
% show_joint(output);
% show_quiver();
% show_gradients();
show_colormaps(stats.u_mean, stats.v_mean, gt, img1, mask, min_image, flowpath);

if exist('flowsubdir','var')
    flowpath = fullfile(flowpath, flowsubdir);
end
save_errors(stats, gt, flowpath);


function show_colormaps(dx, dy, gt, img1, mask, min_image, flowpath)
rows = 1:size(img1,1);
cols = 1:size(img1,2);

alpha = 0.3;
underlay = 255 * double(~mask);
underlay = repmat(alpha * double(underlay(rows, cols)), [1,1,3]);

% figure; clf; colormap(gray(256)); maximize();
%     image(uint8(underlay + dx_overlay)); axis('image');
% figure; clf; maximize();
%     output.var_u(output.var_u > 1.0) = 0;
%     imagesc(log10(output.var_u)); axis('image');
% 
% figure; clf; colormap(gray(256)); maximize();
%     image(uint8(underlay + dy_overlay)); axis('image');
% figure; clf; maximize();
%     output.var_v(output.var_v > 1.0) = 0;
%     imagesc(log10(output.var_v)); axis('image');
% 
% figure; clf; colormap(gray(256)); maximize();
%     image(uint8(underlay + dxy_overlay)); axis('image');

delete(fullfile(flowpath,'*.png'));

% dxy = sqrt(dx.^2 + dy.^2);
% dxy_overlay = (1-alpha) * flow_image(dxy(rows, cols));
% dxy_rgb = 255 * flow_image(dxy(rows, cols), [1 1 1]);

imwrite(uint8(underlay), fullfile(flowpath,'underlay.png'));
imwrite(uint8(min_image), fullfile(flowpath,'min_image.png'));

dx_rgb = 255 * flow_image(dx(rows, cols), [1 1 1]);
imwrite(uint8(dx_rgb), fullfile(flowpath,'u.png'));

dy_rgb = 255 * flow_image(dy(rows, cols), [1 1 1]);
imwrite(uint8(dy_rgb), fullfile(flowpath,'v.png'));

dx_overlay = (1-alpha) * flow_image(dx(rows, cols));
imwrite(uint8(underlay + dx_overlay), fullfile(flowpath,'u_overlay.png'));

dy_overlay = (1-alpha) * flow_image(dy(rows, cols));
imwrite(uint8(underlay + dy_overlay), fullfile(flowpath,'v_overlay.png'));

figure;
    img = [real(gt.flowmap) dx(rows, cols) imag(gt.flowmap) dy(rows, cols)];
    image(flow_image(img)); axis('image');

    
function save_errors(stats, gt, flowpath)

err_u = stats.u_mean - real(gt.flowmap);
mae_u = mean(abs(err_u(~isnan(err_u))));
mean_u_var =       mean(     stats.u_var(~isnan(stats.u_var)));
mean_u_std = real( mean(sqrt(stats.u_var(~isnan(stats.u_var)))) );

err_v = stats.v_mean - imag(gt.flowmap);
mae_v = mean(abs(err_v(~isnan(err_v))));
mean_v_var =       mean(     stats.v_var(~isnan(stats.v_var)));
mean_v_std = real( mean(sqrt(stats.v_var(~isnan(stats.v_var)))) );

disp([mae_u mean_u_var mean_u_std; 
      mae_v mean_v_var mean_v_std]);

datestring = datestr(now,'yyyymmddTHHMMSS');
filename = sprintf('%s_stats.mat', datestring);
parameters = gt.parameters;

save(fullfile(flowpath,filename), ...
     'parameters', 'mae_u', 'mae_v', 'stats');

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
%     fprintf(fid, '\n');
%     fprintf(fid, 'du_mean = %5.3f\n', stats.du_mean);
%     fprintf(fid, 'du_mean_abs = %5.3f\n', stats.du_mean_abs);
%     fprintf(fid, 'du_mean_sqr = %5.3f\n', stats.du_mean_sqr);
%     fprintf(fid, 'du_rms = %5.3f\n', stats.du_rms);
%     fprintf(fid, 'du_var = %5.3f\n', stats.du_var);
%     fprintf(fid, 'du_std = %5.3f\n', stats.du_std);
%     fprintf(fid, '\n');
%     fprintf(fid, 'dv_mean = %5.3f\n', stats.dv_mean);
%     fprintf(fid, 'dv_mean_abs = %5.3f\n', stats.dv_mean_abs);
%     fprintf(fid, 'dv_mean_sqr = %5.3f\n', stats.dv_mean_sqr);
%     fprintf(fid, 'dv_rms = %5.3f\n', stats.dv_rms);
%     fprintf(fid, 'dv_var = %5.3f\n', stats.dv_var);
%     fprintf(fid, 'dv_std = %5.3f\n', stats.dv_std);
    fclose(fid);
end


function show_interior(output, dx, dy)
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

    
function show_joint(output)
n_disp = 101;
max_disp = 0.5;
dx_rng = linspace(-max_disp, max_disp, n_disp);
dy_rng = linspace(-max_disp, max_disp, n_disp);
hough = hist3([output.u(:) output.v(:)], 'edges', {dx_rng, dy_rng});

figure; clf; maximize();
    imagesc(dx_rng, dy_rng, hough);        
    
  
function show_quiver(output, img1, img2)
figure(1); clf; colormap(gray(256)); maximize();
    subplot(2,2,1); imagesc(img1); axis('image');
    subplot(2,2,2); imagesc(img2); axis('image');
    subplot(2,2,3); hold on;
    imagesc(img1);
    quiver(output.x_vec, output.y_vec, output.u_vec, output.v_vec, 5);
    axis('image','ij');
linkaxes([subplot(2,2,1), subplot(2,2,2), subplot(2,2,3)]);        
    

function show_gradients(output)
figure(2); clf; colormap(gray(256)); maximize();
    subplot(2,2,1); imagesc(output.Ix); axis('image');
    subplot(2,2,2); imagesc(output.Iy); axis('image');
    subplot(2,2,3); imagesc(output.It); axis('image');
linkaxes([subplot(2,2,1), subplot(2,2,2), subplot(2,2,3)]);


%% Test script
function test_script()

clc;
% close all;
timebar('closeall');

% imgroot = 'U:\projects\nailfold\capture\2012_10_22\Left\Digit4\x300\';
% % imgpath = fullfile(imgroot, 'seq1\preprocessed\registered_g1d\masked');
% imgpath = fullfile(imgroot, 'seq2\corrected\registered_g1d\masked\halfsize');

imgroot = 'U:\projects\nailfold\synthesis';
d = dir(fullfile(imgroot, '*T*'));
d = d([d(:).isdir]);

% imgroot = fullfile(imgroot, d(end).name);
imgroot = fullfile(imgroot, '20130320T135651');

imgroot = fullfile(imgroot, 'halfsize');

% profile clear; profile on;
pt_flow_fcn(imgroot);
% profile off; profile report;

