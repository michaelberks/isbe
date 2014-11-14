function observationMask = create_observation_mask(imgStackSize, nObservations)
if (nargin==0 && nargout==0), test(); return; end

if (prod(imgStackSize) <= nObservations)
    observationMask = true(imgStackSize);
    return;
end

observationMask = false(imgStackSize);

nPixels = imgStackSize(1)*imgStackSize(2);
nImages = imgStackSize(3);

% Generate map of #observations for each pixel.
% (Initially equally distributed.)
observationMap = floor(nObservations / nPixels) * ones(imgStackSize(1:2));

% Add one to some pixels to get exact number of observations.
nRemaining = nObservations - sum(observationMap(:));
p = randperm(nPixels);
observationMap(p(1:nRemaining)) = observationMap(p(1:nRemaining)) + 1;

for i = 1:imgStackSize(1)
    for j = 1:imgStackSize(2)
        p = randperm(nImages);
        nThisPixel = observationMap(i,j);
        observationMask(i, j, p(1:nThisPixel)) = true;
    end
end


function test()
clc;

imgStackSize = [10,10,10];
nObservations = 0.5 * prod(imgStackSize);

observationMask = create_observation_mask(imgStackSize, nObservations);

% for i = 1:size(observationMask,3)
%     figure(1); clf; hold off; colormap(gray(256));
%         imagesc(observationMask(:,:,i)); axis('image');
% 	pause;
% end

dummy = 0;

