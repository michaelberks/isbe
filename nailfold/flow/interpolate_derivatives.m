function [Ix, Iy, It] = interpolate_derivatives(observations, displacements, ...
                                                f_clear_persistent) 
if (nargin==0 && nargout==0), test(); return; end

if ~exist('f_clear_persistent','var'), f_clear_persistent = false; end

Ix = [];
Iy = [];
It = [];

persistent b;
persistent Ix_padded;
persistent Iy_padded;
persistent It_padded;

if f_clear_persistent, Ix_padded = []; return; end

max_disp = ceil(max(abs(displacements(:))));
if isempty(Ix_padded) || (b < max_disp)
    b = max_disp;
    
    Ix_padded = mb_pad(observations.Ix, [b,b,0], 'replicate');
    Iy_padded = mb_pad(observations.Iy, [b,b,0], 'replicate');
    It_padded = mb_pad(observations.It, [b,b,0], 'replicate');
end

imsz = size(observations.It);
          
Ix = zeros(imsz);
Iy = zeros(imsz);
It = zeros(imsz);

x0 = 1:imsz(2);
y0 = 1:imsz(1);

for i = 1:imsz(3)
    [xx,yy] = meshgrid(x0 - displacements(i,1) + b, ... 
                       y0 - displacements(i,2) + b);
    
    if 0
        % Matlab function with lots of overhead.
        
        Ix(:,:,i) = interp2(Ix_padded(:,:,i), xx, yy, '*linear');
        Iy(:,:,i) = interp2(Iy_padded(:,:,i), xx, yy, '*linear');
        It(:,:,i) = interp2(It_padded(:,:,i), xx, yy, '*linear');
    else
        % My quick and dirty approximation.
        
        Ix(:,:,i) = qinterp2(Ix_padded(:,:,i), xx, yy, '*linear', false);
        % Reuse interpolation points from previous run of qinterp2
        Iy(:,:,i) = qinterp2(Iy_padded(:,:,i), xx, yy, '*linear', true);
        It(:,:,i) = qinterp2(It_padded(:,:,i), xx, yy, '*linear', true);
    end
end


%% Test function
function test()
clc;

imgroot = 'U:\projects\nailfold\synthesis\20130809T170015';

nFrames = 100;
imgStack = load_image_stack(imgroot, nFrames);

[Ix, Iy, It] = image_derivatives(imgStack);
observations.Ix = Ix;
observations.Iy = Iy;
observations.It = It;

gt = load(fullfile(imgroot, '_ground_truth.mat'));
disps = gt.jitter(1:nFrames, :);

interpolate_derivatives([], [], true);
[Ix2, Iy2, It2] = interpolate_derivatives(observations, disps);

for i = 1:nFrames
    figure(1); clf; hold off; colormap(gray(256));
        image(uint8(127+2*Ix2(:,:,i))); axis('image');
    drawnow;
end

return


n = 25;
[xx,yy] = meshgrid(-n:n, -n:n);
rr = sqrt(xx.^2 + yy.^2);

img0 = (rr < 0.7*n) & (rr > 0.3*n);
img0 = 255*double(img0);

f = [1 4 6 4 1];
f = f / sum(abs(f(:)));
img0 = conv2(f,f,img0, 'same');

nF = 10;
j = 1 * randn(nF, 2);
j = conv2(j, [1;2;1], 'same');
j = j - j(ones(nF,1), :);

x = 1:size(img0,2);
y = 1:size(img0,1);

imgStack = zeros([size(img0),nF]);

b = ceil(max(abs(j(:))));
img0 = mb_pad(img0, [b,b,0], 'replicate'); 
for i = 1:nF
    [xx, yy] = meshgrid(x + j(i,1) + b, ...
                        y + j(i,2) + b);
    imgStack(:,:,i) = interp2(img0, xx, yy, '*linear');
end

[Ix, Iy, It] = image_derivatives(imgStack);

observations.Ix = Ix;
observations.Iy = Iy;
observations.It = It;

[Ix2, Iy2, It2] = interpolate_derivatives(observations, j, true);
[Ix2, Iy2, It2] = interpolate_derivatives(observations, j);

for i = 1:nF
    figure(1); clf; hold off; colormap(gray(256));
        image(uint8(127+2*Ix2(:,:,i))); axis('image');
    pause(0.25);
end

