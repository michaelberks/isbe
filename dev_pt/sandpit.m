clc;
clear;
% close all;


warning('off', 'ASYM:unexpectedArgument');
d_args{1}.decomp_type = {'dt'};
d_args{1}.feature_shape = 'rect';
d_args{1}.feature_type = 'complex';
d_args{1}.do_max = 0;
d_args{1}.rotate = 0;
d_args{1}.use_nag = 0;
d_args{1}.levels = 1:3;

d_args{2}.decomp_type = {'gabor'};
d_args{2}.num_angles = 6;
d_args{2}.sigma_range = [1 2 4];    
d_args{2}.do_max = 0;
d_args{2}.rotate = 0;
d_args{2}.feature_type = 'complex';

for i_d = 1:2
    d_args{i_d}.win_size = 3;
    d_args{i_d}.normalise = 0;
    d_args{i_d}.pca = [];
end

a1 = rand(128);%create_gauss_bar(4, 1, 22.5, 128, 128, 64, 64);

for off = 1:8
    a2 = mb_pad(a1, [off off]);

    figure('windowstyle','docked');
    for i_d = 1:2
        r1 = compute_filter_responses(a1, d_args{i_d});
        r2 = compute_filter_responses(a2, d_args{i_d});

        f1 = sample_image_features(r1, 64, 64, d_args{i_d});
        
        i_p = 1;
        for i_x = -1:1
            for i_y = -1:1
                f2 = sample_image_features(r2, 64+off+i_y, 64+off+i_x, d_args{i_d});
                subplot(3,3,i_p); plot(f1, f2, 'x'); hold all;
                axis([-1,1,-1,1],'equal','square');
                i_p = i_p + 1;
            end
        end
    end
    legend({'DT samples', 'Gabor samples'});
end


return

points = dlmread('U:/tmp/ref_pts.txt');
triangles = dlmread('U:/tmp/triangles.txt') + 1;

triangles(49,:) = [51,41,42];
triangles(92,:) = [41,47,42];
triangles(87,:) = [ 3,41,50];
triangles(88,:) = [49, 3,50];

figure(1); clf; hold on;
    for i = 1:size(triangles,1)
        v1 = triangles(i,1);
        v2 = triangles(i,2);
        v3 = triangles(i,3);
        plot(points([v1,v2],1), points([v1,v2],2), 'r-');
        plot(points([v2,v3],1), points([v2,v3],2), 'r-');
        plot(points([v3,v1],1), points([v3,v1],2), 'r-');
        c = mean(points([v1,v2,v3],:));
%         text(c(1), c(2), sprintf('%02d',i-1));
    end
    plot(points(:,1), points(:,2), 'b.');
    b = (max(points(:,1))-min(points(:,1)) ) * 0.1;
    axis('ij', 'equal', ...
         [min(points(:,1))-b, max(points(:,1))+b, ...
          min(points(:,2))-b, max(points(:,2))+b]);
    
f = getframe(gca);

figure(2);
image(f.cdata);
imwrite(f.cdata, 'U:/tmp/face_map.png');

return

imgpath = 'U:\data\nailfolds';
% d = dir([imgpath,'\56409c.png']);
% d = dir([imgpath,'\59214c.png']);
d = dir([imgpath,'\30785c.png']);
axvec = [];
for i = 1
    img = imread(fullfile(imgpath,d(i).name));
    subplot(2,1,1); imshow(img); axvec(end+1) = gca;
    
    % First row is always black (for some unknown reason)
    img = img(2:end,:);
    
    g2d_responses = compute_gaussian_2nd_derivatives(img, [1,2,4]);
    
    scl = 3;
    blobA = 0.5 * (g2d_responses(:,:,scl,2) + g2d_responses(:,:,scl,3));
    blob1 = 0.5 * (g2d_responses(:,:,scl,2) - g2d_responses(:,:,scl,3));
    blob2 = g2d_responses(:,:,scl,1);
    blobB = sqrt(blob1.^2 + blob2.^2);
    
    blob = blobB;
    
    % Try to do away with background
    bg_value = blob(1,1);
    blob(1,:) = nan; blob(:,1) = nan;
    for y = 2:size(blob,1)
        for x = 2:size(blob,2)
            if (blob(y,x) == bg_value) && ...
               (isnan(blob(y-1,x)) || isnan(blob(y,x-1)))
                blob(y,x) = nan;
            end
        end
    end

%     subplot(2,1,2); imagesc(blob); axis('image','off'); axvec(end+1) = gca;
    
    
    a = 0.9;
    img2 = double(repmat(img,[1,1,3]));
    img2 = a*img2;
    
    tA = 0.0; tB = 2;
    mask = (blobA > tA) & (blobB > tB);
    img2(:,:,1) = img2(:,:,1) + (1-a)*double(mask)*255;
    mask = (blobA < -tA) & (blobB > tB);
    img2(:,:,2) = img2(:,:,2) + (1-a)*double(mask)*255;
%     mask = (blobA > 0) && (blobB - abs(blobA) > 0.5);
%     img2(:,:,1) = img2(:,:,1) + (1-a)*double(mask)*255;
    subplot(2,1,2); imshow(img2/255); axvec(end+1) = gca;
    
    figure(); colormap(gray(256));
    subplot(2,1,1); imagesc(blobA);
    subplot(2,1,2); imagesc(blobB);
    
end
linkaxes(axvec);


return

x = linspace(-10,10,201);
vec = [];
smin = 1.0;
smax = 2.0;
for sig = linspace(smin,smax,101)
    cx = cos(2*pi*x/2.35);
%     env = exp(-0.5*(x.^2 / sig^2)); % gaussian
%     env = exp(-abs(x/sig)); % exponential
    env = (1/pi) * (sig ./ (x.^2 + sig*sig)); % cauchy

    f = env .* cx;
    vec(end+1,1:2) = [sig mean(f(:))];

    if (sig == smin || sig == smax)
        figure(); plot(x,f);
    end
end
figure(); plot(vec(:,1),vec(:,2));

return

sigma = 1;

x = linspace(-5*sigma, 5*sigma, 201);
sx = x/sqrt(2*sigma*sigma);
% sy = y/sqrt(2*sigma*sigma);
kernel = -2.205*sx + 0.9780*sx.^3;
gaussian_envelope = exp(-(sx.^2));
f = kernel .* gaussian_envelope;

plot(x,f);

return

t = linspace(-1,1,51);
[xx,yy] = meshgrid(t,t);
zz = xx.*xx.*yy;

figure(1); clf;
mesh(xx,yy,zz);
xyzlabel('x','y','z');


return

x = linspace(-5,5,101);
t = 3;
tb = 2;

y = 0.25*x.^4 + 1.5*x.^2*(t-tb) + 0.75*(t-tb).^2;

figure(1); clf;
plot(x,y);

return

pd = {'a','b','c'};

a = 'd';
switch a
    case pd
        disp('padded');
    otherwise
        disp('unpadded');
end


return

filter_mat = [];
filter_sz = [];
for sigma = [4,2,1]
    filters = gaussian_filters_2d(3,sigma);
    if isempty(filter_sz)
        filter_sz = size(filters);
    end
    
    if (size(filters,1) < filter_sz(1))
        sz_diff = filter_sz(1) - size(filters,1);
        padding = sz_diff/2;
        new_filters = zeros(filter_sz);
        new_filters(1+padding:end-padding, 1+padding:end-padding, :) = ...
            filters;
        filters = new_filters;
    end
    
    filter_mat = [filter_mat reshape(filters,[],3)];
end

for c = 1:size(filter_mat,2)
    filter_mat(:,c) = filter_mat(:,c) - mean(filter_mat(:,c));
end

imgpath = 'U:\projects\mammography\data\retinograms\drive\training\images';
img = imread([imgpath,'\21_training.png']);
img = double(img(:,:,2));

figure(1); clf; colormap(gray(256));
imagesc(img); axis('image','ij');

% [x,y] = ginput(1);
x = 164; y = 417;

patch_width = floor(filter_sz(1) / 2);
patch = img(round(y)-patch_width:round(y)+patch_width,...
            round(x)-patch_width:round(x)+patch_width);
 
figure(2); clf; colormap(gray(256));
imagesc(patch); axis('image','ij');

response0 = filter_mat' * patch(:);
nullspace = null([filter_mat' -response0]);
n_basis = size(nullspace,2);

while true
    coeffs = rand(n_basis,1);
%     coeffs = zeros(n_basis,1); coeffs(ceil(rand*n_basis)) = 1;
    example = nullspace * coeffs;
    example = example / example(end);
    example = reshape(example(1:end-1), filter_sz(1:2));
   
    figure(3); clf; colormap(gray(256));
    imagesc(example); axis('image','ij');
    
    response = filter_mat' * example(:);
    [response0 response]
    
    pause;
end


return

x = linspace(-1,1,201);
y = double(x<0.1);
resp = [];
   
for sigma = 1:40
    [g,dg,ddg,dx] = gaussian_filters_1d(sigma);
    sigma
    
    fy1 = conv2(y,dg,'same');
    fy2 = conv2(y,ddg,'same');
    figure(1); clf; hold on;
        plot(y,'b-');
        plot(fy1,'g-');
        plot(fy2,'r-');
        axis([0,length(x),-0.1,0.1]);
        pause;
    resp = [resp mean(fy(110:111))];
end

figure(2); plot(resp);

return

imlist = create_image_lists({'image_root','A:\asym\data\retinograms\DRIVE',...
                             'image_dir','training'});

return

image_in = ceil(rand(256,256)*255);

inds = ceil(rand(100,1)*numel(image_in));
[rows,cols] = ind2sub(size(image_in), inds);

samples = sample_pixel_data(image_in, rows, cols);

return

[g,dg,ddg] = gaussian_filters_1d(5);
Gx = g'*dg;
Gy = dg'*g;

figure(1); colormap(gray(256));
subplot(1,2,1); imagesc(Gx); axis('image','ij');
subplot(1,2,2); imagesc(Gy); axis('image','ij');

Gxx = Gx(2:end-1,3:end)-Gx(2:end-1,1:end-2);
Gyy = Gy(3:end,2:end-1)-Gy(1:end-2,2:end-1);

G45 = Gy+Gx; G45d = G45(3:end,3:end)-G45(1:end-2,1:end-2);
G135 = Gy-Gx; G135d = G135(3:end,1:end-2)-G135(1:end-2,3:end);

figure(2); colormap(gray(256));
subplot(2,2,1); imagesc(Gxx); axis('image','ij');
subplot(2,2,2); imagesc(Gyy); axis('image','ij');
subplot(2,2,3); imagesc(Gxx-Gyy); axis('image','ij');
subplot(2,2,4); imagesc(G135d-G45d); axis('image','ij');
return

clc; clear;

theta = 60 * pi/180;
halfw = 4;

halfsz = 49;
img = zeros(2*halfsz+1,2*halfsz+1);
v = [-sin(theta); cos(theta)]; % normal to the line
count = 0;
while (count < 100000)
    x = ((2*rand)-1) * halfsz;
    y = ((2*rand)-1) * halfsz;
    if (x*x+y*y > halfsz*halfsz)
        continue;
    end

    p = [x; y];
    d = p'*v;
    if (abs(d) < halfw)
        p = round(p);
        p = halfsz+1 + [p(1), -p(2)];
        img(p(2), p(1)) = img(p(2), p(1)) - 1;
        count = count + 1;
    end
end
img = 2*img / max(abs(img(:))) + 1; % scale to [-1,1]
img = img + 0.1*randn(size(img));
img = 127 + 4*img;

figure(1); clf; hold on; colormap(gray(256));
image(img); axis('image','ij');

sigmas = 2.^linspace(log2(1),log2(8),40);

centre = round(size(img)/2);
centre = [110,82];
centre = [50,50];

responses = zeros(length(sigmas),6);
phi_err = zeros(length(sigmas),1);

for i = 1:length(sigmas)
    sigma = sigmas(i);
    [g,dg,ddg] = gaussian_filters_1d(sigma);

    % define 2D filters in terms of 1D filters
    Gxx = g'*ddg;
    Gxy = dg'*dg;
    Gyy = ddg'*g;

    Ixx = conv2(img, Gxx, 'same');
    Ixy = conv2(img, -Gxy, 'same');
    Iyy = conv2(img, Gyy, 'same');
    
    Ixx = Ixx(centre(1), centre(2));
    Ixy = Ixy(centre(1), centre(2));
    Iyy = Iyy(centre(1), centre(2));
    
    Isum = (Ixx+Iyy)/2;
    Idiff = (Ixx-Iyy)/2;
    
    phi = pi/2 + 0.5 * atan2(2*Ixy, Ixx-Iyy);
    phi_err(i) = abs(ori_error(theta, phi));
%     [theta phi phi_err(i)]
   
    responses(i,:) = [sigma, Ixx, Iyy, Ixy, Isum, Idiff];
end    
disp([responses phi_err]);

figure(2); clf; hold on;
plot(sigmas,responses(:,5) / max(abs(responses(:,5))),'b-');
plot(sigmas, phi_err/max(abs(phi_err)),'r-');
ax = axis;
plot([halfw,halfw],10e6*[-1,1],'k:');
axis(ax);

return

figure(2); colormap('gray');
imagesc(Ixx); axis('image');

x = 24; y = 24;
Ixx = Ixx(y,x); Ixy = Ixy(y,x); Iyy = Iyy(y,x);

t = linspace(0,2*pi,100);
R = Ixx*cos(t).^2 + Iyy*sin(t).^2 + Ixy*sin(2*t);

A = (Ixx+Iyy)/2;
B = sqrt(((Ixx-Iyy)/2)^2 + Ixy*Ixy);
phi = atan2(2*Ixy, Ixx-Iyy);
R2 = A + B*cos(2*t-phi);

figure(3); clf; hold on;
plot(t,R,'b-',t,R2,'r:');

if (A>0),   plot(mod(2*pi + phi/2, pi), 0, 'b*');
else,       plot(mod(2*pi + phi/2+pi/2, pi), 0, 'b*');
end
    
return

halfsz = 25;
thr = halfsz/2;

n = 100000;
A = zeros(n, 2*halfsz+1);
b = zeros(n, 1);
w = zeros(n, 1);
for i = 1:n
    halfw = rand*(halfsz+1);
    w(i) = halfw;
    alpha = rand;
    profile = alpha * tophat(halfsz, halfw) + ...
              (1-alpha) * sine_wave(halfsz, halfw);
    
    A(i,:) = profile;
    
    if (abs(halfw-thr) < 0.5), 
        b(i) = 1;
    else
        b(i) = 0;
    end
end
A = A + 0.02*randn(size(A));
x = A\b;

figure(1); clf; hold on;
rng = -halfsz:halfsz;
plot(rng,x);
ax = axis;
plot(thr*[-1,-1],10e6*[-1,1],'k:');
plot(thr*[1,1],10e6*[-1,1],'k:');
axis(ax);

[g,dg,ddg] = gaussian_filters_1d(thr);
ddg = ddg / max(abs(ddg));
halfsz = (length(ddg)-1)/2;
pad = (length(ddg)-length(x)) / 2;
x = [zeros(pad,1); x; zeros(pad,1)];

n = 100;
halfw = linspace(0, halfsz, n);
responses = zeros(n, 2, 2);
for i = 1:n
    profile = tophat(halfsz, halfw(i));
    responses(i,:,1) = [profile*x -profile*ddg'];

    profile = sine_wave(halfsz, halfw(i));
    responses(i,:,2) = [profile*x -profile*ddg'];
    
    figure(3); clf; hold on;
    plot(-halfsz:halfsz,profile,'b-',...
         -halfsz:halfsz,-ddg,'r-');
    responses(i,2,2)
end
responses = responses(:,:);
for i = 1:size(responses,2)
    responses(:,i) = responses(:,i) / max(abs(responses(:,i)));
end
figure(2); clf; hold on;
plot(halfw, responses);

return

% Approximation to the LinOp filter
theta = 30 * pi/180;
width = 2;

halfw = 7;
img = zeros(2*halfw+1,2*halfw+1);
v = [cos(theta); -sin(theta)];e
for i = 1:1000000
    x = ((2*rand)-1) * halfw;
    y = ((2*rand)-1) * halfw;
    if (x*x+y*y > halfw*halfw)
        continue;
    end

    p = [x; y];
    d = p'*v;
    if (abs(d) < width)
        p = round(p + halfw+1);
        img(p(2), p(1)) = img(p(2), p(1)) + 1;
    end
end

figure(1); clf; hold on; colormap('gray');
imagesc(img); axis('image');

return

clc; 
% clear;

reg = build_vessel_predictor();
real_coeffs = real(reg.beta(2:2:end));
imag_coeffs = imag(reg.beta(2:2:end));

figure(1); clf; hold on;
    plot(real_coeffs,'r.-');
    plot(imag_coeffs,'b.-');
    
return

img = ceil(get(image,'CData')*8);
% close;

figure(1); colormap('gray');
imagesc(img); axis('image');

[g,dg,ddg] = gaussian_filters_1d();

% define 2D filters in terms of 1D filters
Gxx = g'*ddg;
Gxy = dg'*dg;
Gyy = ddg'*g;

Ixx = conv2(img,Gxx,'same');
Ixy = conv2(img,-Gxy,'same');
Iyy = conv2(img,Gyy,'same');

figure(2); colormap('gray');
imagesc(Ixx); axis('image');

x = 24; y = 24;
Ixx = Ixx(y,x); Ixy = Ixy(y,x); Iyy = Iyy(y,x);

t = linspace(0,2*pi,100);
R = Ixx*cos(t).^2 + Iyy*sin(t).^2 + Ixy*sin(2*t);

A = (Ixx+Iyy)/2;
B = sqrt(((Ixx-Iyy)/2)^2 + Ixy*Ixy);
phi = atan2(2*Ixy, Ixx-Iyy);
R2 = A + B*cos(2*t-phi);

figure(3); clf; hold on;
plot(t,R,'b-',t,R2,'r:');

if (A>0),   plot(mod(2*pi + phi/2, pi), 0, 'b*');
else,       plot(mod(2*pi + phi/2+pi/2, pi), 0, 'b*');
end
    

return

imsz = 30;

mask = zeros(2*imsz+1,2*imsz+1);
for i = -imsz:imsz
    for j = -imsz:imsz
        mask(i+imsz+1,j+imsz+1) = (sqrt(i*i+j*j) < imsz-1);
    end
end

linop = -ones(2*imsz+1,2*imsz+1);
w = 3;
linop(imsz-w+1:imsz+w+1,:) = 1;

im2 = -ones(size(linop));

theta = 15 * pi/180;
ct = cos(theta); st = sin(theta);
for i = -imsz:imsz
    for j = -imsz:imsz
        dxy = [ct st; -st ct]*[i;j];
        dx = round(dxy(1)); dy = round(dxy(2));
        if (-imsz < dx) && (dx < imsz) && ...
           (-imsz < dy) && (dy < imsz)
            im2(i+imsz+1,j+imsz+1) = linop(dx+imsz+1,dy+imsz+1);
        end
    end
end

linop = im2;
linop(~mask) = 0;

linop = linop - min(linop(:));
linop = linop / max(linop(:));

figure(1); colormap('gray');
imagesc(linop);
axis('image');

figpath = 'S:\projects\mammography\matlab\papers\9999 journal (orientation)\figs\filtering';
imwrite(linop, [figpath,'/linop.png']);

return

[g,dg,ddg] = gaussian_filters_1d(3);
Gxx = g'*ddg;
Gyy = Gxx';
Gxy = dg'*dg;

figure(1); colormap('gray'); clf;
subplot(1,3,1); imagesc(Gxx+Gyy); axis('image'); 
subplot(1,3,2); imagesc(Gxx-Gyy); axis('image'); 
subplot(1,3,3); imagesc(Gxy); axis('image'); 

disp([rank(Gxx) rank(Gyy) rank(Gxy)]);
disp([rank(Gxx+Gyy) rank(Gxx-Gyy)]);

return

aroot = asymmetryroot;
d_root = [aroot 'data/retinograms/DRIVE/'];

p_root = 'data\models\vessel\width\pc20120229T160934\';
forest = u_load([aroot p_root 'random_forest']);
forest.tree_root = [aroot p_root];

sampling_args = u_load([aroot p_root '\sampling_args.mat']);

sampling_args_c.num_levels = sampling_args.num_levels;
sampling_args_c.feature_shape = sampling_args.feature_shape;
sampling_args_c.feature_type = sampling_args.feature_type;
sampling_args_c.do_max = sampling_args.do_max;
sampling_args_c.rotate = sampling_args.rotate;
sampling_args_c.win_size = sampling_args.win_size;

if isfield(sampling_args, 'use_nag')
    sampling_args_c.use_nag = sampling_args.use_nag;
else
    sampling_args_c.use_nag = args.use_nag;
end

ret_te1 = u_load([d_root 'test/images_extended/01_test_ext.mat']);
f_mask = u_load([d_root 'test/foveal_masks/01_test_f_mask.mat']);

vp = classify_image(...
    'image_in',rgb2gray(ret_te1),... % the mandatory arguments
    'sampling_args',sampling_args_c,...
    'forest', forest, ...
    'decomp_type', 'dt',...
    'forest_type', 'width',...
    'use_probs', 0,...
    'mask', f_mask,...
    'num_trees', [], ...
    'max_size', 128);
figure; imgray(rgb2gray(ret_te1));
figure; imgray(vp);

return

im = imread('U:\projects\nailfold\images\oneloop2.png');
im = double(im(:,:,1));
im = im(2:end,:);

sigma = 4;
if 0
	[g,dg,ddg] = gaussian_filters_1d(sigma);
	Ix = conv2(g,dg,im,'same');
	Iy = conv2(dg,g,im,'same');
	strength = sqrt(Ix.*Ix + Iy.*Iy);
	orientation = atan2(Iy,Ix);
else
	[strength,orientation] = gaussian_clover_line(im,sigma);
end

out = mb_non_maximal_supp(strength,orientation,false,false);

figure(1); clf; colormap(gray(256));
imrgb = repmat(im,[1,1,3]);
imtmp = im; imtmp(out~=0) = 255;
imrgb(:,:,1) = imtmp;
image(uint8(imrgb)); axis('image');


return

retpath = [asymmetryroot('shared'),'data/retinograms/DRIVE/test/'];
load([retpath,'retinogram_properties.mat']);

widths_vec = cat(1,line_widths{:});
max(widths_vec(:))

figure(1); clf;
for i = 1:20
    % load masks
    load([retpath,'vessel_masks\',sprintf('%02d_test_v_mask.mat',i)]);
	
	width_img = zeros(size(vessel_mask));
	width_img(vessel_inds{i}) = line_widths{i};
	subplot(2,1,1); imagesc(width_img>0 & width_img>10); axis('image');

	contrast_img = zeros(size(vessel_mask));
	contrast_img(vessel_inds{i}) = line_contrasts{i};
	subplot(2,1,2); imagesc(vessel_mask + ...
							(contrast_img>0 & contrast_img>60)); axis('image');
end

return


% clc;
clear;

retpath = [asymmetryroot('shared'),'data\retinograms\DRIVE\test\'];
mampath = [asymmetryroot('shared'),'data\synthetic_lines\real512\'];

% default filename format (overwritten where need be)
format_str = '%02d_test_ext_class.mat';

N = 20;
% respath = [retpath,'predictions\g1d\analytic\']; format_str = 'orientations/%03d_test_ori.mat';
% respath = [retpath,'predictions\g1d\linear_regression_1\'];
% respath = [retpath,'predictions\g2d\linear_regression_1\'];
% respath = [retpath,'predictions\mono\linear_regression_1\'];
% respath = [retpath,'predictions\dt\linear_regression_1\'];

N = 100;
respath = [mampath,'predictions\g2d\analytic\']; format_str = 'orientations/%03d_test_ori.mat';

errors_vec = [];
distances_vec = [];
contrasts_vec = [];

for i = 1:20
    % load masks
    load([retpath,'vessel_masks\',sprintf('%02d_test_v_mask.mat',i)]);
    load([retpath,'foveal_masks\',sprintf('%02d_test_f_mask.mat',i)]);
    load([retpath,'vessel_masks\centre_idx\',sprintf('%02d_test_centre_idx.mat',i)]);
    
    % load orientations (ground truth and estimated)
	% note: estimated orientations are scalars (2T) for analytic
	% results, but vectors (cos2T + i.sin2T) for the regressed outputs
    load([retpath,'orientations\',sprintf('%02d_ori1.mat',i)]); % gt_ori
    est_ori = load_uint8([respath,sprintf(format_str,i)]);

    % sort out masks
    vessel_mask(~foveal_mask) = 0;
    centre_mask = vessel_mask;
    centre_mask(centre_mask==1) = vessel_centres;

	% deal with NaNs
    gt_ori(~vessel_mask) = NaN;
    est_ori(~vessel_mask) = NaN;
	vessel_mask = vessel_mask & ~isnan(gt_ori);
    centre_mask = centre_mask & ~isnan(gt_ori);
	
    err_i	= ori_error(est_ori(vessel_mask),gt_ori(vessel_mask));
    err_img = nan(size(vessel_mask));
    err_img(vessel_mask) = err_i;

    % compute line thickness
    img_dt = distance_transform(1e6*vessel_mask);

    % compute line contrast
%     im = imread([retpath,'images_extended\',sprintf('%02d_test.tif',i)]);
%     im = im(:,:,2);
	load([retpath,'images_extended\',sprintf('%02d_test_ext.mat',i)]);
	im = rgb2gray(ret);
    [mn,sd] = local_image_stats(im,7);

    % stack errors, distances and contrasts for centre lines
	used_mask = centre_mask; subset_name = 'centre';
% 	used_mask = vessel_mask; subset_name = 'vessel';
	
	used_inds = find(used_mask);
    errors_vec = [errors_vec; err_img(used_inds)];
    distances_vec = [distances_vec; img_dt(used_inds)];
    contrasts_vec = [contrasts_vec; sd(used_inds)];
end

% debug: compare with errors estiamted by DRIVE_tests.m
load([respath,'errors/orientation_errors.mat']);

errors_vec = abs(errors_vec * 180/pi);
median(errors_vec)
median(abs(prediction_errs) * 180/pi)

figure; clf;
	set(gcf,'name',respath);
subplot(2,2,1);
	img_dt(~used_mask) = 0;
    imagesc(img_dt); 
    axis('image');
subplot(2,2,2);
    sd(~used_mask) = 0;
	imagesc(sd); 
    axis('image');
subplot(2,2,3); cla; hold on;
    plot(distances_vec,errors_vec,'b.');
    unique_dists = unique(distances_vec);
    means_by_dist = zeros(size(unique_dists));
    for i = 1:length(unique_dists)
        sample_inds = (distances_vec==unique_dists(i));
        means_by_dist(i) = mean(errors_vec(sample_inds));
    end
    plot(unique_dists,means_by_dist,'r*');
    axis([0,max(distances_vec)*1.05,0,90]);
subplot(2,2,4); cla; hold on;
    plot(contrasts_vec,errors_vec,'b.');
    [sorted_sd,inds] = sort(contrasts_vec);
    sorted_errs = errors_vec(inds);
    cummean = cumsum(sorted_errs)./(1:length(sorted_errs))';
    plot(sorted_sd,cummean,'r-');
    axis([0,max(contrasts_vec)*1.05,0,90]);
    

figure; clf; hold on;
    plot(distances_vec,errors_vec,'b.','markersize',1);
    unique_dists = unique(distances_vec);
    means_by_dist = zeros(size(unique_dists));
    for i = 1:length(unique_dists)
        sample_inds = (distances_vec==unique_dists(i));
        means_by_dist(i) = mean(errors_vec(sample_inds));
    end
    plot(unique_dists,means_by_dist,'r*');
    axis([0,max(distances_vec)*1.05,0,100]);
	xlabel('Line thickness');
	ylabel('Abs. angular error (degrees)');
	exportfig([respath,'thickness_vs_error-',subset_name]);
figure; clf; hold on;
    plot(contrasts_vec,errors_vec,'b.','markersize',1);
    [sorted_sd,inds] = sort(contrasts_vec);
    sorted_errs = errors_vec(inds);
    cummean = cumsum(sorted_errs)./(1:length(sorted_errs))';
    plot(sorted_sd,cummean,'r-');
    axis([0,max(contrasts_vec)*1.05,0,100]);
	xlabel('Line contrast');
	ylabel('Abs. angular error (degrees)');
	exportfig([respath,'contrast_vs_error-',subset_name]);

return

figure(2); clf;
subplot(1,2,1);
    plot(sd(centre_mask),err_img(centre_mask),'b.');
    axis([0,max(sd(centre_mask))*1.05,0,2]);
subplot(1,2,2);
    plot(sorted_sd,sorted_errs,'b.');
    axis([0,max(sorted_sd)*1.05,0,2]);

return


load('U:\projects\mammography\data\synthetic_lines2\image001.mat');
% test_image = test_image(1:256,1:256);

gaussian_clover_line(test_image,8);

return

err = rand(3,5);
v = (rand(5,1)>0.5);
V = diag(v);

norm(err*V,'fro').^2
v'*err'*err*v
v'*diag(diag(err'*err))*v
trace(err*V*err')

return

clc;

% [xx,yy] = meshgrid(0:7,0:7);
% X = xx.*yy;
% ha = (1:4); hb = ha(end:-1:1);
% ha = [0,1:4,0]; hb = ha(end:-1:1);
% Y = coldfilt(X,ha,hb)
% Y = my_coldfilt(X,ha)
% 
% return

if 0
	img = imread('u:\tmp\pdf_small.jpg');
	img = mean(img,3);
else
	N = 16;
	[xx,yy] = meshgrid(2:N+1,1:N);
	img = mod(xx.*yy,256);
	n_levels = 2;
% 	img = rand(N,N);
end

% f = [0.5; 0; 0.5];
% colfilter(colfilter(img,f)',f)'

% f = [-1; 1];
% my_coldfilt(my_coldfilt(img,f)',f)'

nL = 2;

dt = dtwavexfm2b(img, nL);
my_dt = my_dtwavexfm2b(img, nL);
my_dt{nL}(:,:,1)

return

for i = 1:length(dt)
	err = dt{i}-my_dt{i};
	disp(sum(abs(err(:))));
end

% figure(1); clf; colormap(gray(256));
% 	image(img); axis('image');

return

% load('m:\nailfold\nail_frame.mat');
nail_frame = imread('U:\projects\nailfold\images\n Tonia MooreV1LD4X3LrgMosaic.bmp');

figure(1); colormap(gray);

rng = 101:105;

% [str2,ori2,scl2] = gaussian_2nd_derivative_line(nail_frame,3);
% 	subplot(2,1,2); imagesc(str2); axis('image');
tic;
[str1,ori1,scl1] = gaussian_clover_line(nail_frame,3);
toc;
	subplot(2,1,1); imagesc(str1); axis('image'); 
	
return

str1(rng,rng)
str2(rng,rng)
ori1(rng,rng)-ori2(rng,rng)
scl1(rng,rng)-scl2(rng,rng)

[norm(str1-str2,'fro') norm(ori1-ori2,'fro') norm(scl1-scl2,'fro')]
	
return

%% RF performance on mass_roi
datapath = 'A:\data\mammograms\2004_screening_processed\mass_roi';
fname = '045LML_roi.mat'; xrng = 151:450; yrng = 401:700;
fname = '024RML_roi.mat'; xrng = 151:450; yrng = 401:700;
fname = '039RML_roi.mat'; xrng = 151:450; yrng = 401:700;
fname = '046RML_roi.mat'; xrng = 151:450; yrng = 401:700;
% fname = '060RML_roi.mat'; xrng = 151:450; yrng = 401:700;
% fname = '073LML_roi.mat'; xrng = 151:450; yrng = 401:700;
load(fullfile(datapath,fname));

xrng = 1:size(bg,2); yrng = 1:size(bg,1);

close all;

figure; clf; colormap(gray(256));
	image(bg); axis('image');

bg = bg(yrng,xrng);

figure; clf; colormap(gray(256));
	image(uint8(bg)); axis('image'); colorbar; maximize;

datapath = ['A:\data\orientation_maps\rf\2004_screening_processed\mass_roi'];
orimap = load_uint8(fullfile(datapath,fname));
orimap = orimap(yrng,xrng);
figure; clf; colormap(blueyellow(256));
	image(255*abs(orimap)); axis('image'); colorbar; maximize;
% figure; clf;
% 	image(ori2rgb(orimap)); axis('image'); colorbar; maximize;

[line_orientation, line_map, line_scale, line_strength] = ...
karssemeijer_line_detection(bg,'line_scales',[1,2,4,8]);
scl = max(abs(line_strength(:)));
line_str = (1-line_strength/scl)/2;
figure; clf; colormap(blueyellow(256));
	image(255*line_str); axis('image'); colorbar; maximize;
% figure; clf;
% 	orimap_g2d = line_str.*exp(sqrt(-1)*line_orientation/180*pi);
% 	image(ori2rgb(orimap_g2d)); axis('image'); colorbar; maximize;

return

%% g2d performance on mass_roi
datapath = 'A:\data\mammograms\2004_screening_processed\mass_roi';
load(fullfile(datapath,'045LML_roi.mat'));

% bg = bg(1:200,1:200);
bg = bg(401:700,151:450);

figure(1); clf;
	cmap = [gray(256)];
	image(ind2rgb(uint8(bg),cmap)); axis('image');

% return

scales = [1,2,4,8];
n = length(scales);
max_strength = -inf;
for i = 1:n
	[line_orientation, line_map, line_scale, line_strength] = ...
		karssemeijer_line_detection(bg,...
			'line_scales',[scales(i)]);
		
	strength{i} = -line_strength;
	max_strength = max(max_strength,max(abs(strength{i}(:))));
	
	map{i} = line_map;
end

figure(2); clf;
for i = 1:n
	strmap = 0.5*(1+strength{i}/max_strength);
	cmap = [redgreen(256)];
	subplot(2,2,i); 
		image(ind2rgb(uint8(1+255*strmap),cmap)); axis('off','image');
		title(sprintf('Max = %0.2f',max(abs(strength{i}(:)))));
end

figure(3); clf;
maxmap = map{1}/max(map{1}(:));
for i = 1:n
	maxmap = max(maxmap,map{i}/max(map{i}(:)));
	cmap = [0 0 0; jet(255)];
	subplot(2,2,i);
		image(ind2rgb(uint8(1+255*map{i}/max(map{i}(:))),cmap)); axis('off','image');
		title(sprintf('Max = %0.2f',max(map{i}(:))));
end
figure(4); clf; imagesc(maxmap);

return

button = 0;
figure(2);
while button~=3
	[x,y,button] = ginput(1);
	
	
end

return

%%	look at regions of a mammogram where there are significantly fewer lines
%	pointing to a focal point (i.e. fewer than you would expect for randomly
%	distributed lines).

% method	= 'g2d';
method	= 'rf_thin';

test = '045LML'; p = [508,1672]; 
% 	p = [470,1900]; % rf_thin
% 	p = [779,1489]; % rf_thin
% 	p = [420,1740]; % g2d  

% raw image
inpath	= 'A:\data\mammograms\2004_screening\abnormals\';
load([inpath,test,'.mat']);

% mask
inpath	= 'A:\data\masks\2004_screening\abnormals\';
load([inpath,test,'_mask.mat']);

% orientation map
inpath	= ['A:\data\orientation_maps\',method,'\2004_screening_processed\abnormals\'];

% line map
switch method
	case 'g2d',
		orimap	= load_uint8([inpath,test,'_ori.mat']);
		inpath	= ['A:\data\line_maps\',method,'\2004_screening_processed\abnormals\'];
		load([inpath,test,'_lines.mat']);
		linemap = line_map & mask;
		clear line_map;
	case {'rf','rf_thin'}
		orimap = load_uint8([inpath,test,'_class.mat']);
		linemap = abs(orimap);
% 		inpath	= ['A:\data\line_maps\old_rf\2004_screening_processed\abnormals\'];
% 		linemap = load_uint8([inpath,test,'_class.mat']); 
		orimap = mod(angle(orimap),pi);
		linemap = (linemap>0.5) & mask;
		linemap = linemap .* mask;
end

% compute k-map
spacing = 16; map_sz = 0*spacing;
r_max = 180; r_min = r_max/2; R = r_min/2;

x_min = max(1,p(1)-map_sz); x_max = min(p(1)+map_sz,size(orimap,2));
y_min = max(1,p(2)-map_sz); y_max = min(p(2)+map_sz,size(orimap,1));

figure(2); colormap(gray(256));
	subimg = mammogram(p(2)-r_max:p(2)+r_max,p(1)-r_max:p(1)+r_max);
	image(subimg); axis('off','image');

figure(3); colormap([0 0 0; hsv(255)]);
	oriimg = orimap(p(2)-r_max:p(2)+r_max,p(1)-r_max:p(1)+r_max);
	image(255*oriimg/pi); axis('off','image');

tic;
[f1vec,f2vec,k_mask] = karssemeijer_radial_projection_multiscale(...
					linemap,orimap,...
					'num_angles',1,...
					'x_min',x_min,'x_max',x_max,...
					'y_min',y_min,'y_max',y_max,...
					'r_min',r_min,'r_max',r_max,'R',R,...
					'spacing',spacing,...
					'mask',mask,'debug',true);
toc;
f1	= 0.5 + min(max(f1vec,-50),50)/100;
f2	= 0.5 + min(max(f2vec,-4),4)/8;


return

mappath = 'A:\data\k_stellate_maps\g2d\2004_screening_processed\abnormals\';
mapdir = dir([mappath,'045*ML_f1.mat']);
maskdir = dir([mappath,'045*ML_mask.mat']);

for i = 1
	f1vals = load_uint8([mappath,mapdir(i).name]);
	[f1max,m_inds] = min(f1vals,[],2);
	
	load([mappath,maskdir(i).name]);
	spacing = median(diff(find(mask)));
	
	f1map = nan(size(mask));
	f1map(mask) = f1max;
	f1map = f1map(1:2*spacing:end,1:2*spacing:end);
	f1map(f1map<0) = nan;
	
	indsmap = nan(size(mask));
	indsmap(mask) = m_inds;
	indsmap = indsmap(1:2*spacing:end,1:2*spacing:end);
	indsmap(isnan(f1map)) = nan;
	
	grads(:,:,1) = conv2(f1map,[0,1,-1],'same');
	grads(:,:,2) = conv2(f1map,[-1,1,0],'same');
	grads(:,:,3) = conv2(f1map,[0,1,-1]','same');
	grads(:,:,4) = conv2(f1map,[-1,1,0]','same');
	maxima = all(grads>0,3);
	maxmap = f1map;
	maxmap(~maxima) = nan;
	
	figure(1); clf;
	subplot(2,2,1);	imagesc(indsmap); axis('off','image');
	subplot(2,2,2);	imagesc(f1map); axis('off','image');
	subplot(2,2,4);	imagesc(maxmap); axis('off','image');
end

return

p = 0.5;

Nvec = [10,20,40];

figure(1); clf; 
for iN = 1:length(Nvec)
	N = Nvec(iN);
	k = 1:N;
	for i = 1:length(k)
		pk(i) = nchoosek(N,k(i)) * p^k(i) * (1-p)^(N-k(i));
	end
	mn = p*N;
	sd = sqrt(N*p*(1-p));

	figure(1);
	subplot(length(Nvec),1,iN); hold on;
		plot(k,pk); ax = axis; ax([1,2]) = [0,N];
		plot([mn,mn],[0,1],':');
		plot([mn+sd,mn+sd],[0,1],'--');
		plot([mn-sd,mn-sd],[0,1],'--');
		axis(ax);
end


return

a = mod(1*randn(1,1000),pi);
figure(3);
	hist(a,linspace(0,pi,20));

return

p = rand(1,10); p = p/sum(p(:));
q = rand(1,10); q = q/sum(q(:));

for i = 1:length(p)
	hist_dist(p,q,'match')
	p = p([2:end 1]); q = q([2:end 1]);
end

return

% generate Poisson map

k = 0:50;
kfact = factorial(k);

lambda = 1:max(k);
p_max = zeros(size(lambda));
for i = 1:length(lambda)
	pk = lambda(i).^k ./ kfact .* exp(-lambda(i));
	p_max(i) = max(pk);
end
figure(1); clf;
	plot(lambda,p_max,'b.-');

p_img = p_max'*p_max;
figure(3); clf; hold on;
	imagesc(lambda,lambda,log(p_img)); axis('image','xy');

return

[lp_map,lp_dist,r_map,t_map] = logpolar_map;

figure(1); clf; colormap([0 0 0; jet(max(lp_map(:)))]);
	subplot(2,2,1); image(r_map+1); axis('off','image');
	subplot(2,2,2); image(t_map+1); axis('off','image');
	subplot(2,2,3); image(lp_map+1); axis('off','image');
	
return

N = 10;
a = ceil(N*rand(10000000,1));
% w = rand(size(a));

tic;
h = hist(a,1:N);
toc;

tic;
s = full(sparse(1,a,1, 1,N));
toc;

sum(abs(h-s))

return



f1max = 100;

% test = '002LCC'; mass_c2 = [320,378]; 
% test = '024RCC'; mass_c2 = [400,400]; 
test = '024RML'; mass_c2 = [377,393]; 
% test = '029RCC'; mass_c2 = [400,400]; 
% test = '039RML'; mass_c2 = [460,400]; 
% test = '046RCC'; mass_c2 = [385,340]; 
% test = '056RML'; mass_c2 = [400,400]; 
% test = '060RML'; mass_c2 = [355,425]; mass_c2 = [560,245];
% test = '061RCC'; mass_c2 = [520,535]; % false postive?
% test = '073LML'; mass_c2 = [400,400]; 
% test = '085LCC'; mass_c2 = [400,400]; 
% test = '085LML'; mass_c2 = [400,400]; 
% test = '097RML'; mass_c2 = [400,380]; 
% test = '101RCC'; mass_c2 = [400,400]; 
% test = '112LCC'; mass_c2 = [400,440]; 

method = 'g2d';

metapath = 'A:\data\mammograms\2004_screening\abnormals\meta\';
metadir = dir([metapath,'*_meta.mat']);

for i = 1:length(metadir)
	test = metadir(i).name(1:6)
	
	f1path = ['A:\data\k_stellate_maps\',method,'\2004_screening_processed\abnormals\'];
	f1 = load_uint8([f1path,test,'_f1.mat']);
	f1 = 0.5 + min(max(f1,-f1max),f1max)/(2*f1max); % [-f1max,f1max] -> [0,1]

	load([f1path,test,'_mask.mat']);
	mask = mask(1:2:end,1:2:end);

	map = nan(size(mask));
	map(mask) = f1(:,3);

	load([metapath,test,'_meta.mat'],'meta_xy');
	meta_xy = meta_xy/2;
	meta_xy(end+1,:) = meta_xy(1,:);

	cmap = [0 0 0; redgreen(255)];
	cmap = jet(256);
	figure(1); clf; hold on; colormap(cmap);
		image(uint8(1+255*map)); axis('image','off','ij');
		plot(meta_xy(:,1),meta_xy(:,2),'k-');

	outpath = ['U:\projects\mammography\data\k_stellate_maps\',method,'\'];
	exportfig([outpath,test,'_f1'],{'eps','pdf','png'});
end	


return


inpath = 'M:\asymmetry_project\Paper submissions\2011 BMVC (Orientation)\matlab_figs\';
outpath = 'S:\projects\mammography\matlab\papers\2011bmvc\figs\';

figdir = dir([inpath,'*.fig']);
for i = 1:length(figdir)
	open([inpath,figdir(i).name]);
% 	graph; 
	exportfig([outpath,figdir(i).name]);
end

return

% imhw	= 3;
% imsz	= 2*imhw+1;
% img		= reshape(1:imsz*imsz,[imsz,imsz])

imhw	= 5;
[xx,yy] = meshgrid(-imhw:imhw,-imhw:imhw);
img		= yy;

iimg00	= integral_image(img);
iimg45	= integral_image_diag(img);
iimg45	= integral_image_diag(img')';

x = imhw+1; y = x;
inds_str = haar_inds('diag',3);
haar1 = gen_haar_image(inds_str);
inds_str = haar_inds('diag45',3);
haar2 = gen_haar_image(inds_str);

figure(1); colormap(gray(256));
subplot(1,2,1);	imagesc(haar1); axis('image');
subplot(1,2,2);	imagesc(haar2); axis('image');



return

% normal square feature
w = 1; h = w;
% inds	= [	h w 1; 0 w -1; h 0 -1; 0 0 1];
inds	= [	0 -1 1; -h -h-1 -1; w -w-1 -1; w-h -w-h-1 1;
			0 0 1; -w w -1; h h -1; h-w w+h 1];
haar = gen_haar_image(inds,'diag_pos');
inds	= [	-1 0 1; -h-1 -h -1; -w-1 w -1; -w-h-1 w-h 1;
			0 0 1; w -w -1; h h -1; w+h h-w 1];
haar1 = haar - gen_haar_image(inds,'diag_neg');

% normal square feature
w = round(w*sqrt(2)); h = w; 
haar2 = gen_haar_image(inds,'normal');
border = (size(haar1,1)-size(haar2,1))/2;
haar2 = [zeros(border,size(haar2,2)); haar2; zeros(border,size(haar2,2))];
haar2 = [zeros(size(haar2,1),border) haar2 zeros(size(haar2,1),border)];

figure(1); colormap(gray(256));
subplot(1,2,1);	imagesc(haar1); axis('image');
subplot(1,2,2); imagesc(haar2); axis('image');

return

% % get sum over square of width 4 at centre of image (should be 16)
% x = 3; y = 3; w = 4; h = 4;
% sum =	  iimg00(y+h-1,x+w-1) ...
% 		- iimg00(y+h-1,x-1) ...
% 		- iimg00(y-1,  x+w-1) ...
% 		+ iimg00(y-1,  x-1)
% 
% % get sum over square of width 4 at centre of image (should be 25)
% iimg45	= integral_image_diag(img)
% x = 2; y = 4; w = 4; h = 4;
% sum =	  iimg45(y+w,   x+w) ...
% 		+ iimg45(y+h,   x-h) ...
% 		- iimg45(y,     x) ...
% 		- iimg45(y+w+h, x+w-h)

return

% test logistic regression
overlap = 0.2;
X = 2*rand(1000,2);
y = (X(:,1)>1);
X(y,1) = X(y,1)-overlap/2;
X(~y,1) = X(~y,1) + overlap/2;
t = rand*2*pi;
X = X * [cos(t) -sin(t); sin(t) cos(t)]';

logreg = pt_log_reg_train(X,y);
y_fit = logistic_regressor_predict(logreg,X);

cmap = hsv(512); cmap = cmap(1:256,:);
figure(1); clf; colormap(cmap); hold on;
subplot(2,1,1);
	scatter(X(:,1),X(:,2),16,y); axis('equal');
subplot(2,1,2); 
	scatter(X(:,1),X(:,2),16,y_fit); axis('equal');
	
	
return

x = linspace(-3,3,101);
[xx,yy] = meshgrid(x,x);
r = sqrt(xx.*xx + yy.*yy);

G = exp(-0.5*r.*r);
figure(1); clf; colormap(gray(256));
	imagesc(G); axis('image');
	
dG = -r.*G;
figure(2); clf; colormap(gray(256));
	imagesc(dG); axis('image');


return

x = linspace(-3,3,101);
[xx,yy] = meshgrid(x,x);

r = sqrt(xx.*xx+yy.*yy);
s1 = 1; s2 = 2;
H1 = (s1^-2)*exp(-0.5*r.*r*s1^-2);
H2 = (s2^-2)*exp(-0.5*r.*r*s2^-2);
H = H1 - H2;

figure(1); clf; sbsz = [1,3];
mysubplot(sbsz); imagesc(H1);
mysubplot(sbsz); imagesc(H2);
mysubplot(sbsz); imagesc(H);

figure(2); clf;
mesh(H);

return

x	= linspace(0,2*pi,101);
y	= sin(x) + 0.1*randn(size(x));

% filters
fL	= exp(-0.5*(-3:3).^2); 
fL	= fL/sum(fL);
fH	= zeros(1,7); fH(4) = 1;
fH	= fH-fL;

yL = conv(y,fL); yL = yL(4:end-3);
yH = conv(y,fH); yH = yH(4:end-3);
figure(1); clf; hold on;
	plot(x,yL,'r-',x,yH,'g-',x,yL+yH,'b-');

display(norm(y-yL-yH,'fro'));

return

figure(1); clf; hold on;
k = 0:50;
mn = 4;

lm = mn;
p = lm.^k*exp(-lm)./factorial(k);
plot(k,p,'b.-');

sig = sqrt(mn);
p = 1/sqrt(2*pi*sig*sig)*exp(-0.5 * (k-mn).^2/(sig*sig));
plot(k,p,'r.-');


return

clear tree;

% load most recently created forest if none specified
forest_root		= [asymmetryroot,'data/line_orientation_rfs/'];
forest_dir		= dir([forest_root,'pc*']);
forest_job		= [forest_dir(end).name,'/']
forest_dir		= dir([forest_root,forest_job,'/random_forest*.mat']);
forest_fname	= [forest_root,forest_job,forest_dir(end).name];
forest			= u_load(forest_fname);

% get sampling arguments
sampling_args	= u_load([forest_root,forest_job,'/sampling_args.mat']);

% generate some test data using same sampling arguments
fprintf('Generating test data...');
sampling_args = rmfield(sampling_args,'detection_type');
sampling_args.num_samples = 8000;
sampling_args.task_id = 1;
[testX,testy] = sample_saved_dt_line_data(sampling_args);
fprintf('done\n');

fprintf('Predicting...');

% predict with linear regressor
for t = 1:length(forest.trees)
	load([forest.tree_root,forest.tree_dir,sprintf('traindata%04d.mat',t)]);
	R = angle(y)\X;
	y_fit = exp(sqrt(-1)*testX*R');
	e_linear(:,t) = ori_error(testy,y_fit) * 180/pi;
end

% predict with tree
y_fit = mb_random_forest_reg_predict(forest,testX);
% y_fit = abs(y_fit).*exp(sqrt(-1)*2*angle(y_fit));
e_tree = ori_error(testy,y_fit) * 180/pi;

fprintf('done\n');

errvec = mean(e_tree,2);
	err_stats.mean = mean(errvec);
	err_stats.sd = std(errvec);
errvec = sort(abs(errvec));
	err_stats.abs_mean = mean(errvec);
	err_stats.abs_median = median(errvec);
	err_stats.abs_range = errvec([1,end])';
	pcs = round([0.01:0.01:1]*length(errvec));
	err_stats.abs_percentiles = errvec(pcs)';
display(err_stats);

errvec = mean(e_linear,2);
	err_stats.mean = mean(errvec);
	err_stats.sd = std(errvec);
errvec = sort(abs(errvec));
	err_stats.abs_mean = mean(errvec);
	err_stats.abs_median = median(errvec);
	err_stats.abs_range = errvec([1,end])';
	pcs = round([0.01:0.01:1]*length(errvec));
	err_stats.abs_percentiles = errvec(pcs)';
display(err_stats);

return

fprintf('...');
fprintf('done\n');

for t = 1:3
	newtree = mb_tree_reg_train(X,y,...
		'random_m',forest.d,...
		'split_min',10);
	figure; show_tree_output(newtree);
end

return
newtree = tree;

leaves = find(newtree.var==0);
for L = 1:length(leaves)
	figure(1); show_tree_output(newtree);
	hold on;
	plot(newtree.outputs{leaves(L)},'b.');
	newtree.nodeerr(leaves(L))
	abs(mean(newtree.outputs{leaves(L)}))
	pause;
end

leaves = find(newtree.var==0);
disps = 1-newtree.nodeerr(leaves);
figure(2); hist(disps,0.4:0.01:1);

return

load('U:\projects\mammography\data\line_orientation_rfs\pc20110218T171730\01_trees\traindata.mat');

tree_args.random_m = forest.d;
tree_args.split_min = 50;
tree_args.prune = 0;
tree_args.names = {};
tree_args.mod = 0;
tree_args.w_prior = 0.01;

rand('twister',1234);
newtree = mb_tree_reg_train(X, y, tree_args);

tmp = [newtree.var newtree.goodness]; [tmp(1:3,:) tmp(1:3,2)./tmp(1:3,3)]

leaves = find(newtree.var==0);
for L = 1:length(leaves)
	figure(1); show_tree_output(newtree);
	hold on;
	plot(newtree.outputs{leaves(L)},'b.');
	newtree.nodeerr(leaves(L))
	abs(mean(newtree.outputs{leaves(L)}))
	pause;
end

leaves = find(newtree.var==0);
disps = 1-newtree.nodeerr(leaves);
figure(2); hist(disps,0.4:0.01:1);


return


N = 1000;
t = linspace(0,2*pi,N)';
X = [cos(t) sin(t)];

meanvec = zeros(length(t),2);
for i = 1:length(t)-1
	D = X(1:i,:); 
	res = sum(D,1)/size(D,1); 
	meanvec(i,1) = sqrt(res*res');
	D = X(i+1:end,:); 
	res = sum(D,1)/size(D,1); 
	meanvec(i,2) = sqrt(res*res');
end
d = meanvec(:,2)-meanvec(:,1);
d = -d.*d;

figure(1); clf; hold on;
	plot(	t,meanvec(:,1),'b-',...
			t,meanvec(:,2),'r-',...
			t,d,'g-');
	

return


X(1,:) = sort(rand(1,N));
X(2,:) = rand(1,N);
z = zeros(1,N);
y = zeros(1,N);

t = 1*linspace(-1,1,101);
p = -tanh(t)*pi;

inds = ceil(X(1,:)*length(t));
y = complex( cos(p(inds)),sin(p(inds)) );

zz = angle(y);

figure(10); clf; hold on; colormap(redblue);
scatter(X(1,:),X(2,:),24,zz,'.');
axis([0,1,0,1],'square');

tree_args = struct(	'random_m', 0,...
										'split_min', 10, ...
										'prune', 0,...
										'mod', 0, ...
										'do_circular', [],...
										'impure_thresh', 1e-6,...
										'names', []);
tree = mb_tree_reg_train(X', y', tree_args);

return

clear forest;

% load([asymmetryroot,'data\line_orientation_rfs/pc20110211T143503/random_forest01.mat']);
% forest = random_forest;

if ~exist('forest','var')
	forest_root		= [asymmetryroot,'data/line_orientation_rfs/'];
	
	forest_dir		= dir([forest_root,'pc*']);
	forest_job		= [forest_dir(end).name,'/'];
	
	forest_dir		= dir([forest_root,forest_job,'/random_forest*.mat']);
	forest_fname	= [forest_root,forest_job,forest_dir(end).name];
	forest				= u_load(forest_fname);
	
	% load training data
	load([forest_root,forest_job,'01_trees/traindata.mat']);
end

in		= X(:,1:30).*exp(complex(0,X(:,31:end)));
in		= reshape(in,[size(in,1),6,5]);
out		= angle(y);

% set up labels for variables
labels = {};
for level = 1:5
	for subband = 1:6
		labels{subband,level} = sprintf('L%i:S%i',level,subband);
	end
end

ncols	= 256; t = linspace(0,2*pi,ncols);
cmap	= [(cos(t)+1)/2; (sin(t)+1)/2; zeros(1,ncols)]';

% default parameters
axlim = [];

%% get inputs in pairs
pair = 'ph_mag';
switch lower(pair)
	case 're_im'
		% real vs imag
		in = in(:,:); labels = labels(1:end);
		in = cat(3,real(in),imag(in)); 
		labels = cat(1,labels,labels);
		prefix = cat(1,repmat({'Re'},[1,size(in,2)]),repmat({'Im'},[1,size(in,2)]));
		fname = 'dt_scatter_Re_Im';
		sbsz = [5,6];
		
	case 'ph_mag',
		% phase vs magnitude
% 		in = in(:,:); labels = labels(1:end);
% 		in = cat(3,angle(in),log(abs(in)));
		in = in(:,:); labels = labels(1:end);
		in1 = angle(in); in2 = log(abs(in));
		labels = cat(1,labels,labels);
		prefix = cat(1,repmat({'Ph'},[1,size(in,2)]),repmat({'Mag'},[1,size(in,2)]));
		fname = 'dt_scatter_Ph_Mag';
		sbsz = [5,6];
		axlim = [-pi,pi,-3,3];

	case 'mag_mag',
		% mag vs mag, phase vs phase
		in = in(:,[1,6,2,5,3,4],:); in = in(:,:); % put similar orientations next to each other
		labels = labels([1,6,2,5,3,4],:); labels = labels(1:end);
		in = cat(3,[abs(in(:,1:2:end)) angle(in(:,1:2:end))],...
							 [abs(in(:,2:2:end)) angle(in(:,2:2:end))]); 
		in = reshape(in,[size(in,1)*3,size(in,2)/3,2]);
		in = [in(:,1:size(in,2)/2,:); in(:,1+size(in,2)/2:end,:)];
		in = reshape(in,[size(in,1)/6,size(in,2)*6,2]);
		labels = cat(2,reshape(labels,[2,length(labels)/2]),reshape(labels,[2,length(labels)/2]));
		prefix = cat(1,repmat({'Mag'},[6,size(in,2)/6]),repmat({'Ph'},[6,size(in,2)/6]));
		prefix = reshape(prefix,[2,size(in,2)]);
		fname = 'dt_scatter_Mag_Mag';
		sbsz = [5,6];

	case 'mag_mag2',
		% mag vs mag, phase vs phase
		in1 = [abs(in(:,:,2:end)) 2*angle(in(:,:,2:end))]; labels1 = labels(:,2:end);
		in2 = [abs(in(:,:,1:end-1)) angle(in(:,:,1:end-1))]; labels2 = labels(:,1:end-1);
		
		sbsz = [8,6];
	
		fname = 'dt_scatter_Mag_Mag2';
		
	case 'relmag',
		% mag vs mag, phase vs phase
		in1 = [log(abs(in(:,[2,3,4,5,6,1],:))./abs(in(:,[1,2,3,4,5,6],:)))]; labels1 = labels(:,2:end);
		in2 = [angle(in)]; labels2 = labels(:,1:end-1);
		
		sbsz = [5,6];
		axlim = [-5,5,-pi,pi];
	
		fname = 'dt_scatter_RelMag';
end

figure(1); clf; colormap(cmap);
d = 1; inds = reshape(1:prod(sbsz),sbsz(end:-1:1))';
for x1 = 1:sbsz(1)
	for x2 = 1:sbsz(2)
		mysubplot(sbsz,inds(x1,x2)); hold on;
		scatter(in1(:,d), in2(:,d), 4, out, '.');
		if ~isempty(axlim), axis(axlim); end
% 		axis('equal'); 
		set(gca,'box','on');
% 		plot(1e6*[-1,1],1e6*[-1,1],'-','color',0.7*[1,1,1]);
% 		xlabel(sprintf('%s(%s)',prefix{1,d},labels{1,d}));
% 		ylabel(sprintf('%s(%s)',prefix{2,d},labels{2,d}));
		d = d+1;
	end
end

set(1,'paperposition',[0,0,40,30]);
exportfig([figpath,fname],{'fig','eps','pdf','png'});

return


D = forest.D;
varhst = zeros(100,D);

vardisp = zeros(100,D);
vardispcount = zeros(100,D);

parmat = zeros(D,D);

width = [];
contrast = [];
orientation = [];
squash = [];

max_level = 0;

% figure(1); clf; hold on;
for t = 1:length(forest.trees)
	% load parameters
	load([forest.tree_root,forest.tree_dir,...
				sprintf('../line_parameters/01/parameters%03d.mat',t)]);
	width = [width [parameters.width]];
	contrast = [contrast [parameters.contrast]];
	orientation = [orientation [parameters.orientation]];
	squash = [squash [parameters.squash]];
	
			
	% load tree
	load([forest.tree_root,forest.tree_dir,forest.trees{t}]);
	leaves = find(all(tree.children==0,2));

	% sort leaves in order of output
	[sorted,inds] = sort(angle(tree.class(leaves))/2);
	leaves = leaves(inds);

% 	figure(1);
% 		for L = 1:length(leaves)
% 			leaf = leaves(L);
% 			n = length(tree.outputs{leaf});
% 			plot(L*ones(1,n),angle(tree.outputs{leaf})/2,'r.');
% 		end
% 		plot(angle(tree.class(leaves))/2);
% 	axis([ 0,length(leaves)+1,pi/2*[-1,1] ]);
% 	exportfig('U:\matlab\figs\mammography\scatter');

% 	varhst = varhst + hist(tree.var(tree.var~=0),60);

	% compute level of the tree for each node
	tree.level = nan(size(tree.node));
	tree.level(1) = 1;
	for i = 2:length(tree.level)
		tree.level(i) = tree.level(tree.parent(i))+1;
	end
	for i = 1:length(tree.level)
		if tree.var(i)>0
			% branch node - count
			varhst(tree.level(i),tree.var(i)) = ...
				varhst(tree.level(i),tree.var(i)) + 1;
		else
			% leaf node - compute dispersal
			
			% trace all the way back to root, adding dispersal value to every
			% ancestor of node i
			par = tree.parent(i);
			while par>0
				vardisp(tree.level(par),tree.var(par)) = ...
					vardisp(tree.level(par),tree.var(par)) + 1-tree.nodeerr(i);
				vardispcount(tree.level(par),tree.var(par)) = ...
					vardispcount(tree.level(par),tree.var(par)) + 1;
				par = tree.parent(par);
			end
		end
	end
	max_level = max(max_level,max(tree.level)-1);
	
	% fill in co-occurence matrix for neighbouring nodes in tree
	for node = 2:length(tree.node)
		v1 = tree.var(node); % cut variable 1
		if (v1==0), continue; end
		v2 = tree.var(tree.parent(node)); % cut variable 2
		parmat(v1,v2) = parmat(v1,v2) + 1;
	end
end

imagesc(parmat+parmat'); colormap(gray(256)); axis('image');
exportfig([figpath,'neighbours']);

return

varhst = varhst(1:max_level,:);

% normalize dispersal values
vardisp = vardisp./vardispcount;
vardisp = vardisp(1:max_level,:);

figure(2); clf; colormap(gray);
	sbsz = [3,1];
	im1 = varhst./(max(varhst,[],2)*ones(1,D));
	im2 = varhst./(ones(size(varhst,1),1)*max(varhst,[],1));
% 	im1 = vardisp./(max(vardisp,[],2)*ones(1,D));
% 	im2 = vardisp./(ones(max_level,1)*max(vardisp,[],1));
	mysubplot(sbsz,1); imagesc(im1); xlim([0,D+1]);
	mysubplot(sbsz,2); imagesc(im2); xlim([0,D+1]);
	mysubplot(sbsz,3); hold on;
		bar(1:D,sum(varhst,1));
% 		plot(1:D,sum(varhst,1));
		plot([1;1]*(0:6:D)+0.5,[0;1e6]*ones(1,length(0:6:D)),'k:');
		axis([0,D+1,0,max(sum(varhst,1))]);
		

figure; set(gcf,'paperposition',[0,0,20,10]);
	colormap(gray(256));
	clf; 
		imagesc(im1); xlim([0,D+1]); 
		xlabel('Feature'); ylabel('Tree level');
		exportfig([figpath,'histim1']);
	clf; 
		imagesc(im2); xlim([0,D+1]); 
		xlabel('Feature'); ylabel('Tree level');
		exportfig([figpath,'histim2']);
	clf; hold on;
		colormap('default');
		bar(1:D,sum(varhst,1));
% 		plot(1:D,sum(varhst,1));
		plot([1;1]*(0:6:D)+0.5,[0;1e6]*ones(1,length(0:6:D)),'k:');
		axis([0,D+1,0,max(sum(varhst,1))]);
		xlabel('Feature'); ylabel('Count');
		exportfig([figpath,'feathist']);
	clf; hold off;
		h3 = reshape(sum(varhst,1),[6,10]);
		bar3(h3);
		xlabel('Feature'); ylabel('Orientation subband');
		exportfig([figpath,'bar3d']);
	close('all');

return

sbsz = [2,2]; nh = 24;
figure(3); clf;
	mysubplot(sbsz); hist(contrast,nh); title('contrast');
	mysubplot(sbsz); hist(width,nh); title('width');
	mysubplot(sbsz); hist(orientation,nh); title('ori');
	mysubplot(sbsz); hist(squash,nh); title('squash');


    