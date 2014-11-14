function pyramid = build_image_pyramid(imgStack, nLevels, smoother)
% Build an image pyramid of the imageStack. The output is a cell array with
% an imageStack, of successively smaller dimensions, in each cell.

if (nargin==0 && nargout==0), test(); return; end

if ~exist('nLevels','var') || isempty(nLevels), nLevels = 1; end

% The kernel used to smooth the image in each dimension.
if ~exist('smoother','var') || isempty(smoother), smoother = [1,2,1]/4; end

% Pad image stack so that we don't lose pixels at lower levels.
multiple = 2^(nLevels-1);
m = size(imgStack, 1);
n = size(imgStack, 2);
imsz = [m,n];
pad_size = ceil(imsz/multiple)*multiple - imsz;
imgStack = mb_pad(imgStack, [pad_size, 0], 'replicate', 'post');

% Warn if non-unit smoother.
if (abs(sum(smoother)-1) > 1e-12)
    warning('Smoother does not sum to unity');
end

% First level is always the original image stack.
pyramid{1} = imgStack;

% Return if nLevels is 1.
if (nLevels == 1)
    return
end

for i = 2:nLevels
    pad_size = zeros(1,3);
    for d = 1:2
        pad_size(d) = (length(smoother)-1) / 2; % for smoothing
    end
    imgStack = mb_pad(imgStack, pad_size, 'replicate');

    for d = 1:2
        % Ignore identity smoother
        if (smoother == 1)
            continue;
        end

        % Reshape the smoother.
        dims = ones(1,3);
        dims(d) = length(smoother);
        smoother = reshape(smoother, dims);

        imgStack = convn(imgStack, smoother, 'valid');
    end
    
    imgStack = imgStack(1:2:end, 1:2:end, :);
    
    pyramid{i} = imgStack;
end


%% Test function
function test()
clc;

% Get Easter egg image that comes with image().
h = figure; image();
img = get(get(gca,'children'),'cdata');
close(h);

img = img(1:57, 1:57);
imgStack = repmat(img, [1,1,3]);

nLevels = 4;
pyramid = build_image_pyramid(imgStack, nLevels);

figure(1); clf; colormap(gray(256));
for i = 1:nLevels
    subplot(1,nLevels,i);
    imagesc(uint8(pyramid{i}(:,:,1)));
    axis('image');
end


