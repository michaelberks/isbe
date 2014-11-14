% PROBABLY JUNK NOW

function output = pt_flow_derivs_joint(Ix, Iy, It, ...
                                       step, lambda, tau)
% Use image_derivatives() to get Ix, Iy and It.

if (nargin==0 && nargout==0), test_script(); return; end

if ~exist('step','var') || isempty(step), step = 1; end
if ~exist('lambda','var') || isempty(lambda), lambda = 0.0; end % Regularizer
if ~exist('tau','var') || isempty(tau), tau = 0.0; end

n_pixels = numel(It);
imginds = reshape(1:n_pixels, size(It));

[ny,nx] = size(It);
imrng_x = 1:step:nx;
imrng_y = 1:step:ny;

nPatches = length(imrng_x) * length(imrng_y);

nConsData = nPatches;
bigUdata = spalloc(nConsData, n_pixels, nConsData);
bigVdata = spalloc(nConsData, n_pixels, nConsData);
bigTdata = zeros(nConsData, 1);
dataRows = 1;

% Make these persistent as they don't change from one frame to the next
% (unless you use a data-dependent penalty, e.g. steerable random field)
persistent imsz;
persistent bigUsmooth;
persistent bigVsmooth;
persistent bigTsmooth;

% if isempty(imsz) || ...
%    numel(imsz) ~= numel(size(It)) || ...
%    any(imsz ~= size(It))
    imsz = size(It);    
    [bigUsmooth, bigVsmooth, bigTsmooth] = ...
        create_spatial(imsz, imrng_x, imrng_y, tau/6, tau/12);
% end

%% Create the data matrices
for x = imrng_x
    for y = imrng_y
        c = imginds(y, x);
        
        % Data constraints
        bigUdata(dataRows, c) = Ix(y, x);
        bigVdata(dataRows, c) = Iy(y, x);
        bigTdata(dataRows, 1) = It(y, x);
        
        dataRows = dataRows + length(dataRows);
    end
end

z = sparse(size(bigUsmooth,1), size(bigUsmooth,2));
A = [bigUdata       bigVdata; 
     bigUsmooth     z;
     z              bigVsmooth];
b = [bigTdata; 
     bigTsmooth;
     bigTsmooth];

AtA = A'*A + lambda*speye(2*n_pixels, 2*n_pixels);
Atb = A'*b;

% uv = A \ b;
uv = AtA \ Atb;
% uv = pcg(AtA, Atb, [], 100);

uv = reshape(uv, [n_pixels, 2]);

u = full(reshape(uv(:,1), size(It)));
v = full(reshape(uv(:,2), size(It)));

rss = nan(size(It));
var_u = nan(size(It));
var_v = nan(size(It));
interior = nan(size(It));

output = struct(... % Inputs and workspace variables
                'Ix', Ix, 'Iy', Iy, 'It', It, ...
                'x_rng', imrng_x, 'y_rng', imrng_y, ...
                'interior', interior, ...
                ... % Outputs
                'u', u, 'v', v, ...
                ... % Confidence
                'rss', rss, 'var_u', var_u, 'var_v', var_v, ...
                'Us', bigUsmooth, 'Vs', bigVsmooth, 'Ts', bigTsmooth);


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

It = img0(2:end,   2:end);
img2 = img0(1:end-1, 1:end-1);

patch_hw = 0;
step = 1;
lambda = 0.01;
tau = 1.0;

profile clear; profile on;
output = pt_flow_horn_schunk(It, img2, [], [], patch_hw, step, lambda, tau);
profile report; profile off;

% wu = 1 ./ output.var_u;
wu = 1;
u = wu .* output.u; 

% wv = 1 ./ output.var_v;
wv = 1;
v = wv .* output.v;

figure(1); clf; colormap(gray(256));
    subplot(3,2,1); imagesc(It); axis('image','ij');
    subplot(3,2,2); imagesc(img2); axis('image','ij');
    subplot(3,2,3); image(flow_image(u)); axis('image','ij');
    subplot(3,2,4); image(flow_image(v)); axis('image','ij');
    subplot(3,2,5); imagesc(output.var_u); axis('image','ij');
    subplot(3,2,6); imagesc(output.var_v); axis('image','ij');
    
 disp([sum(u(:))/sum(wu(:)) sum(v(:))/sum(wv(:))]);
