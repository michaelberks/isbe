% PROBABLY JUNK NOW

function output = pt_flow_derivs(img1, img2, Ix1, Iy1, ...
                                 patch_hw, step, lambda)

if (nargin==0 && nargout==0), test_script(); return; end

if ~exist('Ix1', 'var'), Ix1 = []; end
if ~exist('Iy1', 'var'), Iy1 = []; end
if ~exist('lambda','var'), lambda = 0.0; end % Regularizer
if ~exist('step','var'), step = 1; end
if ~exist('patch_hw','var'), patch_hw = 3; end

% f = isbe_fspecial('gaussian', [3,3], 0.5);
% img1 = conv2(img1, f, 'same');
% img2 = conv2(img2, f, 'same');

% Use precomputed values of Ix and Iy (e.g. from a previous iteration).
if isempty(Ix1)
    Ix1 = conv2(img1, 0.5*[1,0,-1], 'same');
    Ix1(:,[1,end]) = 0; % border effects
end
if isempty(Iy1)
    Iy1 = conv2(img1, 0.5*[1,0,-1]', 'same');
    Iy1([1,end],:) = 0; % border effects
end

Ix2 = conv2(img2, 0.5*[1,0,-1], 'same');
Ix2(:,[1,end]) = 0; % border effects

Iy2 = conv2(img2, 0.5*[1,0,-1]', 'same');
Iy2([1,end],:) = 0; % border effects

It = (img2 - img1);

imrng_y = (patch_hw+1) : step : (size(img1,1)-patch_hw);
imrng_x = (patch_hw+1) : step : (size(img1,2)-patch_hw);

% Regression parameters
% Pixel weights
if 1
    w = ones(2*patch_hw+1);
else
    w = isbe_fspecial('gaussian', (2*patch_hw)+1 * [1,1], (patch_hw+0.5)/2);
end
w = w / sum(abs(w(:)));
Wmat0 = diag(w(:));

% Workspace variables
imsz = size(img1);
dx_img = nan(imsz);
dy_img = nan(imsz);
var_u = nan(imsz);
var_v = nan(imsz);
rss = nan(imsz);
interior = false(imsz);

pixel_index = 0;
for x = imrng_x
    xrng = (x-patch_hw):(x+patch_hw);
    
    for y = imrng_y
        pixel_index = pixel_index + 1;
        
        % Get average flow over square patch
        yrng = (y-patch_hw):(y+patch_hw);
        
        patch = img1(yrng, xrng);
        if any(any(patch==0)), continue; end
        
        px = Ix1(yrng, xrng);
        py = Iy1(yrng, xrng);

        Wmat = Wmat0;
        if 0
            grey_diff = (patch - patch(patch_hw+1, patch_hw+1));
            Wgrey = exp(-8 * grey_diff.*grey_diff);
            Wgrey = Wgrey / sum(abs(Wgrey(:)));
            Wmat = Wmat .* diag(Wgrey(:));
        end
        
        % Lucas-Kanade solution:
        %                  [px(:) py(:)] * [u v] = [-pt(:)]
        % [px(:) py(:)]' * [px(:) py(:)] * [u v] = [px(:) py(:)]' * [-pt(:)]
        X = [px(:) py(:)];

%         XtX = X' * (Wmat * X);
        XtX = X' * X;
        
        % Check condition number of matrix
        % condXtX = cond(XtX); % slow
        d = XtX(1)*XtX(4) - XtX(3)*XtX(2);
        if (d == 0), continue; end
        
        s = svd(XtX); 
        condXtX = s(1)/s(2); % faster
        if (condXtX > 1e8), continue; end
        
        pt = It(yrng, xrng);
        t = -pt(:);
%         XtY = X' * (Wmat * t);
        XtY = X' * t;

        beta = (XtX + lambda*eye(2)) \ (XtY);
        
        dx_img(y, x) = beta(1);
        dy_img(y, x) = beta(2);
        
        y_hat = X * beta;
        residual = y_hat - t;
        rss(y, x) = residual' * residual;
        
        [N, p] = size(X);
        sigma_hat_sqr = rss(y, x) / (N - p - 1);
        var_beta = XtX \ [sigma_hat_sqr 0; 0 sigma_hat_sqr];
        var_u(y, x) = var_beta(1, 1);
        var_v(y, x) = var_beta(2, 2);

        interior(y, x) = true;
    end
end

output = struct(... % Inputs and workspace variables
                'img', img1, ...
                'Ix', Ix2, 'Iy', Iy2, 'It', It, ...
                'x_rng', imrng_x, 'y_rng', imrng_y, ...
                'interior', interior, ...
                ... % Outputs
                'u', dx_img, 'v', dy_img, ...
                ... % Confidence
                'rss', rss, 'var_u', var_u, 'var_v', var_v);

            
%% Test script            
function test_script()
clc;

hw = 31; % half width
fw = 2*hw + 1; % full width

patch_hw = 5;

img0 = zeros(fw,fw);
% img0(:, hw:end) = 1; % horizontal edge
% img0(hw:end, :) = 1; % vertical edge
img0(hw:end, hw:end) = 1; % corner

% % Cone
% for i = 0:hw
%     for j = 0:hw
%         r = sqrt(i*i + j*j);
% %         value = -r;
%         value = 255*double(r<(hw/2));
%         
%         img0(hw+i+1, hw+j+1) = value;
%         img0(hw-i+1, hw+j+1) = value;
%         img0(hw+i+1, hw-j+1) = value;
%         img0(hw-i+1, hw-j+1) = value;
%     end
% end
% % img0 = img0-1;

if 1
    % f = isbe_fspecial('gaussian', [9,9], 3);
    s = 5;
    f = ones(s,s); f = f / sum(abs(f(:)));
    img0 = conv2(img0, f, 'valid');
    img0 = img0 + 0.05*randn(size(img0));
end

img1 = img0(2:end,   2:end);
img2 = img0(1:end-1, 1:end-1);

patch_hw = 1;
step = 1;
lambda = 0.0;
output = pt_flow_derivs(img1, img2, [], [], patch_hw, step, lambda);

cmap = [zeros(1,3); redgreen(254)];

% wu = 1 ./ output.var_u;
wu = 1;
u = wu .* output.u; 
u = 1+ceil(254 * normim(u(patch_hw+2:end-patch_hw-1, patch_hw+2:end-patch_hw-1), 'stretch_fixed'));
u_rgb = ind2rgb(uint8(u), cmap);

% wv = 1 ./ output.var_v;
wv = 1;
v = wv .* output.v;
v = 1+ceil(254 * normim(v(patch_hw+2:end-patch_hw-1, patch_hw+2:end-patch_hw-1), 'stretch_fixed'));
v_rgb = ind2rgb(uint8(v), cmap);

figure(1); clf; colormap(gray(256));
    subplot(3,2,1); imagesc(img1); axis('image','ij');
    subplot(3,2,2); imagesc(img2); axis('image','ij');
    subplot(3,2,3); imagesc(uint8(255*u_rgb)); axis('image','ij');
    subplot(3,2,4); imagesc(uint8(255*v_rgb)); axis('image','ij');
    subplot(3,2,5); imagesc(output.var_u); axis('image','ij');
    subplot(3,2,6); imagesc(output.var_v); axis('image','ij');
    
 disp([sum(u(:))/sum(wu(:)) sum(v(:))/sum(wv(:))]);
