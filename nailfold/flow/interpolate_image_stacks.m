function stacks_out = interpolate_image_stacks(stacks_in, ...
                                               displacements, ...
                                               frame_indices, ...
                                               f_clear_persistent)
% Compute interpolated image stacks, where each image in every stack is
% shifted by the corresponding row in <displacements>.
% Corresponding images in each of the stacks_in are displaced by the same
% vector.
                                           
if (nargin==0 && nargout==0), test(); return; end

if ~exist('frame_indices', 'var'), frame_indices = 1:size(stacks_in{1},3); end
if ~exist('f_clear_persistent','var'), f_clear_persistent = false; end

persistent b;
persistent padded_stacks;

if f_clear_persistent
    clear('b','padded_stacks');
    return; 
end
if isempty(stacks_in), return; end

max_disp = ceil(max(abs(displacements(:))));
if isempty(padded_stacks) || (b < max_disp)
    b = max_disp;

    padded_stacks = cell(size(stacks_in));
    for s = 1:numel(stacks_in)
        padded_stacks{s} = mb_pad(stacks_in{s}, [b,b,0], 'replicate');
    end
end

% Initialize stacks_out with stacks_in - will be overwritten later
for s = 1:numel(stacks_in)
    stacks_out{s} = stacks_in{s}(:,:,frame_indices);
end

[ymax, xmax, nImages] = size(stacks_in{1});
x0 = 1:xmax;
y0 = 1:ymax;

for i = 1:length(frame_indices)
    [xx,yy] = meshgrid(x0 - displacements(frame_indices(i),1) + b, ... 
                       y0 - displacements(frame_indices(i),2) + b);
    
    if 0
        % Matlab function with lots of overhead.
        for s = 1:numel(stacks_in)
            stacks_out{s}(:,:,i) = interp2(padded_stacks{s}(:,:,frame_indices(i)), ...
                                           xx, yy, '*linear');
        end
        
    else
        % My quick and dirty approximation.
        f_recycle = false;
        for s = 1:numel(stacks_in)
            stacks_out{s}(:,:,i) = qinterp2(padded_stacks{s}(:,:,frame_indices(i)), ...
                                            xx, yy, '*linear', f_recycle);
                                        
            % Reuse workspace variables for subsequent interpolations.
            f_recycle = true;
        end
    end
end


%% Test function
function test()
clc;

imgroot = 'U:\projects\nailfold\synthesis\20130809T170109';

nFrames = 100;
imgStack = load_image_stack(imgroot, nFrames);

[Ix, Iy, It] = image_derivatives(imgStack);

gt = load(fullfile(imgroot, '_ground_truth.mat'));
disps = gt.jitter(1:nFrames, :);
disps = zeros(size(disps));

interpolate_image_stacks([], [], true);
stacks_out = interpolate_image_stacks({imgStack, Ix}, disps);

for i = 1:nFrames
    figure(1); clf; hold off; colormap(gray(256));
        imagesc(stacks_out{2}(:,:,i)); axis('image');
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

