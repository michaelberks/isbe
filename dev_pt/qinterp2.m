function zi = qinterp2(z, xx, yy, method, f_recycle)
% Quick approximation to interp2()

if (nargin==0 && nargout==0), test(); return; end

if ~exist('method','var'), method = ''; end
if ~exist('f_recycle','var'), f_recycle = false; end

xi = xx(1,:);
yi = yy(:,1);
zi = z;

persistent xf;
persistent xc;
persistent xf_wts;
persistent xc_wts;

if ~f_recycle
    % Recompute indices and weights
    xf = floor(xi);
    xc = ceil(xi);
    f_wts = xc - xi;
    c_wts = 1 - f_wts;
    
    [m,ignore] = size(zi);
    xf_wts = f_wts(ones(m,1),:);
    xc_wts = c_wts(ones(m,1),:);
end

% Interpolate z along x
zi = zi(:,xf) .* xf_wts + ...
     zi(:,xc) .* xc_wts;
 
persistent yf;
persistent yc;
persistent yf_wts;
persistent yc_wts;

if ~f_recycle
    % Recompute indices and weights
    yf = floor(yi);
    yc = ceil(yi);
    f_wts = yc - yi;
    c_wts = 1 - f_wts;
    
    [ignore,n] = size(zi);
    yf_wts = f_wts(:,ones(1,n));
    yc_wts = c_wts(:,ones(1,n));
end
    
% Interpolate zi along x
zi = zi(yf,:) .* yf_wts + ...
     zi(yc,:) .* yc_wts;


%% Test function
function test()
clc;

% Get Easter egg image that comes with image().
h = figure; image();
img = get(get(gca,'children'),'cdata');
close(h);

xi = 1:0.25:size(img,2);
yi = 1:0.25:size(img,1);
[xx,yy] = meshgrid(xi, yi);

N = 200;

profile clear; profile on;
tic; for i = 1:N, zi = interp2(img, xx, yy); end; toc;
tic; for i = 1:N, zi2 = qinterp2(img, xx, yy); end; toc;
tic; for i = 1:N, zi2 = qinterp2(img, xx, yy, [], true); end; toc;
profile off; profile report;

err = max(abs(zi(:)-zi2(:)));
disp(err);

figure(1); clf; colormap(gray(256));
    subplot(2,2,1); imagesc(img); axis('image');
    subplot(2,2,2); imagesc(img); axis('image');
    subplot(2,2,3); imagesc(zi); axis('image');
    subplot(2,2,4); imagesc(zi2); axis('image');
    