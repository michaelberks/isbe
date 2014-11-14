% PROBABLY JUNK NOW

clc;

imgroot = 'U:\projects\nailfold\synthesis\';
d = dir(fullfile(imgroot, '*T*'));
imgroot = fullfile(imgroot, d(end).name);

gt = load(fullfile(imgroot, '_ground_truth.mat'));

%% Find sampling positions
xstep = gt.xx(1,2) - gt.xx(1,1);
ystep = gt.yy(2,1) - gt.yy(1,1);

% Remember that xx and yy have been reversed
valid_pts = (gt.pts(:,1) > gt.xx(end,end)) & ...
            (gt.pts(:,1) < gt.xx(1,1)) & ...
            (gt.pts(:,2) > gt.yy(end,end)) & ...
            (gt.pts(:,2) < gt.yy(1,1));
        
pt_inds = find(valid_pts);
n_pts = length(pt_inds);

xinds = zeros(1,n_pts);
yinds = zeros(1,n_pts);

for i = 1:n_pts
    ind = pt_inds(i);
    
    dx = abs(gt.xx(1,:) - gt.pts(ind,1));
%     xinds = find(dx < xstep);
%     xwts(i) = 1 - dx(xinds)/xstep;
% 
    dy = abs(gt.yy(:,1) - gt.pts(ind,2));
%     yinds = find(dy < ystep);
%     ywts = 1 - dy(yinds)/ystep;

    [ignore, xinds(i)] = min(dx);
    [ignore, yinds(i)] = min(dy);
end

%% Load the images and sample along the centreline
d = dir(fullfile(imgroot, 'frame_*.png'));
% d = d(1:200);

n_frames = length(d);

vals = zeros(n_pts, length(d));
for i = 1:n_frames
    img = imread(fullfile(imgroot, d(i).name));
    
    if ~exist('img_inds','var')
        img_inds = sub2ind(size(img), yinds, xinds);
    end
    vals(:,i) = img(img_inds);
end

%% Normalize the sampled points
avg = mean(vals, 2);
for i = 1:n_frames
    mc = [vals(:,i) ones(n_pts,1)] \ avg;
    vals(:,i) = mc(1)*vals(:,i) + mc(2);
    vals(:,i) = vals(:,i) ./ avg;
end
vals0 = vals - mean(vals(:));

f_rng = 1:n_frames;
% f_rng = 1:150;
% f_rng = 171:320;
vals = vals0(f_rng,:);

figure(1); clf; colormap(gray(256));
    imagesc(vals); axis('image');

[g,dg,ddg] = gaussian_filters_1d(2);
Ixx =  conv2(g', ddg, vals, 'same');
Ixy = -conv2(dg', dg, vals, 'same');
Iyy =  conv2(ddg', g, vals, 'same');

response = Ixx+Iyy;
response(response > 0) = nan;
response = abs(response);

figure(2); clf; colormap(gray());
    imagesc(response); axis('image');

% theta_im = 0.5 * atan2(Ixx-Iyy, 2*Ixy);
% theta_im = 0.5 * atan(2*Ixy ./ (Ixx-Iyy));
theta_im = 0.5 * atan2(2*Ixy, Ixx-Iyy);
inds =  (~isnan(response)) & ...
        (response > max(response(:))*0.25);
theta_vec = theta_im(inds);

inds =  (~isnan(response));
wts = response(inds);
wts = wts / sum(abs(wts));

theta_est = sum(wts .* theta_im(inds));
disp(theta_est * [1 180/pi]);

theta_im(isnan(response)) = NaN;
theta_im = theta_im;% + pi/2;
theta_cpx = complex(cos(2*theta_im), sin(2*theta_im));
figure(4);
    image(complex2rgb(theta_cpx)); axis('image');
    
figure(3); clf;
    hist(theta_vec,101);
    
return

figure(1); hold on;
    xc = size(vals,2) / 2;
    yc = size(vals,1) / 2;
    L = 100;
    t = [theta_est] + pi/2;
    for i = 1:length(t)
        plot(xc + [-1,1]*L*cos(t(i)), yc + [-1,1]*L*sin(t(i)), 'r-');
    end

[x,y] = ginput(2);
dx = x(2)-x(1);
dy = y(2)-y(1);
t_tagged = atan2(dy, dx);
disp((t_tagged-pi/2) * [1 180/pi]);
