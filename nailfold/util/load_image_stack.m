function imgStack = load_image_stack(imgroot, t_range, x_range, y_range)
% Load an image sequence and store them in a 3D matrix

if (nargin==0 && nargout==0), test(); return; end

if ~exist('t_range', 'var'), t_range = []; end
if ~exist('x_range', 'var'), x_range = []; end
if ~exist('y_range', 'var'), y_range = []; end

d = dir(fullfile(imgroot, 'frame_*.png'));
if isempty(d)
    error(['Directory ', imgroot, ' does not exist or is empty.']);
    return; 
end

% Compute range for a scalar value
if (isempty(t_range))
    t_range = 1:length(d);
elseif (numel(t_range) == 1)
    t_range = 1:t_range;
end

% Select only the desired frames.
if (0 < min(t_range(:))) && (max(t_range(:)) <= length(d))
    d = d(t_range);
end

% Get the size of each frame
img = mean(imread(fullfile(imgroot, d(1).name)), 3);

if isempty(x_range)
    x_range = 1:size(img, 2);
end
if isempty(y_range)
    y_range = 1:size(img, 1);
end

imgStack = zeros([numel(y_range), numel(x_range), numel(t_range)]);

for i = 1:length(d)
    img = mean(imread(fullfile(imgroot, d(i).name)), 3);
    imgStack(:,:,i) = img(y_range, x_range);
end


%% Test function
function test()
clc;
imgpath = 'U:\projects\nailfold\capture\2012_10_16\Left.Digit4.x300\seq2';

imgStack = load_image_stack(imgpath);
disp(size(imgStack));

imgStack = load_image_stack(imgpath, 11:20);
disp(size(imgStack));

imgStack = load_image_stack(imgpath, [], 161:480, 121:360);
disp(size(imgStack));

for i = 1:size(imgStack,3)
    figure(1); clf; colormap(gray(256));
        image(uint8(imgStack(:,:,i)));
        axis('image');
    pause(0.05);
end
