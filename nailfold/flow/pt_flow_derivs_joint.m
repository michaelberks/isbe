% PROBABLY JUNK NOW

function output = pt_flow_derivs_joint(img1, img2, Ix1, Iy1, ...
                                 patch_hw, step, lambda, tau)

if (nargin==0 && nargout==0), test_script(); return; end

if ~exist('Ix1', 'var'), Ix1 = []; end
if ~exist('Iy1', 'var'), Iy1 = []; end
if ~exist('patch_hw','var'), patch_hw = 0; end
if ~exist('step','var') || isempty(step), step = 1; end
if ~exist('lambda','var') || isempty(lambda), lambda = 0.0; end % Regularizer
if ~exist('tau','var') || isempty(tau), tau = 0.0; end

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

It = (img1 - img2); % = -(img2 - img1)

n_pixels = numel(img1);
imginds = reshape(1:n_pixels, size(img1));

ny = size(img1,1);
imrng_y = max(2,patch_hw+1) : step : min(ny-1,ny-patch_hw);
nx = size(img1,2);
imrng_x = max(2,patch_hw+1) : step : min(nx-1,nx-patch_hw);

nPatches = length(imrng_x) * length(imrng_y);

nConsData = nPatches * (2*patch_hw + 1)^2;
bigUdata = spalloc(nConsData, n_pixels, nConsData);
bigVdata = spalloc(nConsData, n_pixels, nConsData);
bigTdata = zeros(nConsData, 1);
dataRows = 1:(2*patch_hw + 1)^2;

persistent imsz;
persistent bigUspatial;
persistent bigVspatial;
persistent bigTspatial;

if isempty(imsz) || any(imsz ~= size(img1))
    imsz = size(img1);    
    [bigUspatial, bigVspatial, bigTspatial] = ...
        create_spatial(imsz, imrng_x, imrng_y, tau);
end

%% Create the data matrices
for x = imrng_x
    xrng = (x-patch_hw):(x+patch_hw);
    
    for y = imrng_y
        % Get average flow over square patch
        yrng = (y-patch_hw):(y+patch_hw);
        
        dx = Ix1(yrng, xrng);
        dy = Iy1(yrng, xrng);
        dt = It(yrng, xrng);
        
        c = imginds(y, x);
        
        % Data constraints
        bigUdata(dataRows, c) = dx(:);
        bigVdata(dataRows, c) = dy(:);
        bigTdata(dataRows, 1) = dt(:);
        
        dataRows = dataRows + length(dataRows);
    end
end

A = [bigUdata    bigVdata; 
     bigUspatial bigVspatial];
b = [bigTdata; 
     bigTspatial];

AtA = A'*A + lambda*speye(2*n_pixels, 2*n_pixels);
Atb = A'*b;

% uv = A \ b;
uv = AtA \ Atb;
% uv = pcg(AtA, Atb, [], 100);

uv = reshape(uv, [n_pixels, 2]);

u = full(reshape(uv(:,1), size(img1)));
v = full(reshape(uv(:,2), size(img1)));

rss = nan(size(img1));
var_u = nan(size(img1));
var_v = nan(size(img1));
interior = nan(size(img1));

output = struct(... % Inputs and workspace variables
                'img', img1, ...
                'Ix', Ix2, 'Iy', Iy2, 'It', It, ...
                'x_rng', imrng_x, 'y_rng', imrng_y, ...
                'interior', interior, ...
                ... % Outputs
                'u', u, 'v', v, ...
                ... % Confidence
                'rss', rss, 'var_u', var_u, 'var_v', var_v, ...
                'Us', bigUspatial, 'Vs', bigVspatial, 'Ts', bigTspatial);


function [Us, Vs, Ts] = create_spatial(imsz, imrng_x, imrng_y, tau)

if (tau == 0)
    Us = [];
    Vs = [];
    Ts = [];
    return;
end

nPixels = prod(imsz);
imginds = reshape(1:nPixels, imsz);

nPatches = length(imrng_x) * length(imrng_y);

nConsSpatial = nPatches * 4;
Us = spalloc(nConsSpatial, nPixels, 2*nConsSpatial);
Vs = spalloc(nConsSpatial, nPixels, 2*nConsSpatial);
Ts = spalloc(nConsSpatial, 1, 0);

spatialRows = 1:4;

for x = imrng_x
    for y = imrng_y
        % Get average flow over square patch
        n_inds = imginds(y-1:y+1, x-1:x+1);
        
        Us(spatialRows, n_inds(2,2)) = tau;
        Vs(spatialRows, n_inds(2,2)) = tau;

        spatialCols = [             n_inds(1,2) ...
                        n_inds(2,1)             n_inds(2,3) ...
                                    n_inds(3,2)             ];
        inds = (spatialCols-1)*size(Us,1) + spatialRows;
        Us(inds) = -tau;
        Vs(inds) = -tau;
        
        spatialRows = spatialRows + length(spatialRows);
    end
end


%% Test script            
function test_script()
clc;

hw = 71; % half width
fw = 2*hw + 1; % full width

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

img0 = rand(175,73);

if 1
    sigma_n = 0.05;
    
    % f = isbe_fspecial('gaussian', [9,9], 3);
    s = 3;
    f = ones(s,s); f = f / sum(abs(f(:)));
    img0 = conv2(img0, f, 'valid');
    img0 = img0 + sigma_n*randn(size(img0));
end

img1 = img0(2:end,   2:end);
img2 = img0(1:end-1, 1:end-1);

patch_hw = 0;
step = 1;
lambda = 0.01;
tau = 1.0;

profile clear; profile on;
output = pt_flow_derivs_joint(img1, img2, [], [], patch_hw, step, lambda, tau);
profile report; profile off;

% wu = 1 ./ output.var_u;
wu = 1;
u = wu .* output.u; 

% wv = 1 ./ output.var_v;
wv = 1;
v = wv .* output.v;

figure(1); clf; colormap(gray(256));
    subplot(3,2,1); imagesc(img1); axis('image','ij');
    subplot(3,2,2); imagesc(img2); axis('image','ij');
    subplot(3,2,3); image(flow_image(u)); axis('image','ij');
    subplot(3,2,4); image(flow_image(v)); axis('image','ij');
    subplot(3,2,5); imagesc(output.var_u); axis('image','ij');
    subplot(3,2,6); imagesc(output.var_v); axis('image','ij');
    
 disp([sum(u(:))/sum(wu(:)) sum(v(:))/sum(wv(:))]);
