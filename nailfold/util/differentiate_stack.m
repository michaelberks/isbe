function [Ix, Iy, It] = differentiate_stack(imgStack, imgStackNext)
% Return image derivatives in x, y and t (time) for a given image stack.

if (nargin==0 && nargout==0), test(); return; end

if ~exist('imgStackNext','var'), imgStackNext = []; end

smoothers{1} = [1 2 1];
smoothers{2} = [1 2 1];
smoothers{3} = [1];

pad_size = zeros(1,3);
for d = 1:3
    pad_size(d) = (length(smoothers{d})-1) / 2; % for smoothing
    pad_size(d) = pad_size(d) + 1; % for differentiation
end

% Initialize smoothedStack with padded imgStack.
% Superfluous values will be removed during filtering.
smoothedStack = mb_pad(imgStack, pad_size, 'replicate');
if ~isempty(imgStackNext)
    smoothedStackNext = mb_pad(imgStackNext, pad_size, 'replicate');
end    

for d = 1:3
    smoother = smoothers{d};
    
    if (numel(smoother) == 1)
        % Identity smoother
        continue;
    end
    
    smoother = smoother / sum(abs(smoother(:)));

    dims = ones(1,3);
    dims(d) = length(smoother);
    smoother = reshape(smoother, dims);
    
    smoothedStack = convn(smoothedStack, smoother, 'valid');
    if ~isempty(imgStackNext)
        smoothedStackNext = convn(smoothedStackNext, smoother, 'valid');
    end
end

% Differentiate.
filter = [1,0,-1] / 2;
Ix = convn(smoothedStack, reshape(filter, [1,3,1]), 'valid');
Iy = convn(smoothedStack, reshape(filter, [3,1,1]), 'valid');

if isempty(imgStackNext)
    % Compute temporal differences from consecutive frames.
    It = convn(smoothedStack, reshape(filter, [1,1,3]), 'valid');

    % Could also be written as:
%     It = convn(smoothedStack, reshape([1,-1], [1,1,2]), 'valid');
%     It = convn(It, reshape([1,1] / 2, [1,1,2]), 'valid');
else
    % Compute temporal differences between smoothedStack and smoothedStackNext.
    % (This is used in hierarchical methods where the temporal difference
    % is between an image and a warped version of it.)
    It = smoothedStackNext(:,:,2:end) - smoothedStack(:,:,1:end-1);
    
    % Average temporal derivative at a frame and the frame before it to
    % make this measurement comparable to that computed from consecutive
    % frames.
    It = convn(It, reshape([1,1]/2, [1,1,2]), 'valid');
end

% Trim the superfluous borders.
Ix = mb_pad(Ix, [-1, 0, -1]);
Iy = mb_pad(Iy, [0, -1, -1]);
It = mb_pad(It, [-1, -1, 0]);

% Compensate for underestimation at borders
Ix(:,[1,end],:) = 2 * Ix(:,[1,end],:);
Iy([1,end],:,:) = 2 * Iy([1,end],:,:);
It(:,:,[1,end]) = 2 * It(:,:,[1,end]);


%% Test script
function test()
clc;

M = rand(3,3,10);
[Mx, My, Mz] = differentiate_stack(M);

[Mx, My, Mz2] = differentiate_stack(M, M); % Should give the same result

Mz(:,:,3)
Mz2(:,:,3)


