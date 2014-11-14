function imgStackOut = warp_pyramid(imgStack, flowmap)

if (nargin==0 && nargout==0), test(); return; end

uu = real(flowmap);
vv = imag(flowmap);

uu(isnan(uu)) = 0;
vv(isnan(vv)) = 0;

% Allocate space for interpolated imgStack.
imgStackOut = zeros(size(imgStack));

% Compute points at which to sample images
[xx,yy] = meshgrid(1:size(imgStack,2), ...
                   1:size(imgStack,1));
xx = xx + uu;
yy = yy + vv;

% Pad the border with replicated values
border = ceil(max([abs(uu(:)); abs(vv(:))]));
imgStack = mb_pad(imgStack, [border,border,0], 'replicate');
xx = xx + border;
yy = yy + border;

% Interpolate every frame (probably rather slow)
for f = 1:size(imgStack,3)
    imgStackOut(:,:,f) = interp2(imgStack(:,:,f), xx, yy, 'linear');
end


function test()
clc;

imgpath = 'U:\projects\nailfold\synthesis\20130906T133337';
imgStack = load_image_stack(imgpath, 1:10,[],[]);

gt = load(fullfile(imgpath, '_ground_truth.mat'));
flowmap = gt.flowStackMean;

warpedStack = warp_pyramid(imgStack, flowmap);

figure(1); clf; colormap(gray(256));
    imagesc([imgStack(:,:,1) warpedStack(:,:,2)]);
    axis('image','off');






