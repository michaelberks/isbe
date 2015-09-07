function [flowPyramidEst, flowConfidence] = ...
    estimate_flow_multilevel(imageRoot, nFrames, outroot, levels)
% Estimate the flow field in an image sequence using a hierarchical
% estimation framework. Flow is estimated at the coarsest level, then the
% flow field is used to 'warp' the images before computing the image
% gradients at the next level, and so on.

if (nargin==0 && nargout==0), test(); return; end

if isnumeric(imageRoot)
    imgStack = imageRoot;
    nFrames = size(imgStack, 3);
else
    imgStack = load_image_stack(imageRoot, nFrames);
end

nPyramidLevels = levels(end);
imgPyramid = build_image_pyramid(imgStack, nPyramidLevels);

% Start at coarsest level of the pyramid
level = nPyramidLevels;

% Estimate flow at top level (lowest resolution) of the pyramid.
% p_opt is the optimized parameter vector (ie the flow field)
% p_init is the initial parameter vector
% confidence is a confidence map based on the Hessian at every pixel
gt = load_ground_truth(imageRoot, nPyramidLevels, level, nFrames);
outputPath = fullfile(outroot, sprintf('level_%02d', level));
[p_opt, p_init, confidence] = ...
    estimate_flow(imgPyramid{level}, outputPath, gt, []);

flowPyramidEst{level} = complex(p_opt.uu, p_opt.vv);
flowConfidence{level} = confidence.uu .* confidence.vv;

for level = levels(end-1:-1:1)
    % Upsample estimated flow.
    % (Flow magnitude is also doubled at this scale.)
    p_est = [];
    p_est.flow = upscale_flow(flowPyramidEst{level+1}, 1);
    
    % Upscale estimated displacements from previous level
    if isfield(p_opt, 'displacements') && ~isempty(p_opt.displacements)
        p_est.displacements = 2 * p_opt.displacements;
    end
    
    % Estimate update to flow
    % p_est is the estimated flow at the current level.
    gt = load_ground_truth(imageRoot, nPyramidLevels, level, nFrames);
    outputPath = fullfile(outroot, sprintf('level_%02d', level));
    [p_opt, p_init, confidence] = ...
        estimate_flow(imgPyramid{level}, outputPath, gt, p_est);

    flowPyramidEst{level} = p_est.flow + complex(p_opt.uu, p_opt.vv);
    flowConfidence{level} = confidence.uu .* confidence.vv;
end


function out = upscale_flow(in, ntimes)

% ntimes = 0 => no change
if (ntimes == 0)
    out = in;
    return
end

scl = 2^ntimes;
if 0
    out = kron(in, scl*ones(scl,scl));
else
    in = scl * mb_pad(in, 1, 'replicate', 'post');
    out = interp2(in, ntimes);
    out = out(1:end-1, 1:end-1, :);
end


%% Test function
function test()
clc;
clear variables;

% Ensure that the working directory is that which contains the modified
% stdnls.m and color.m files so that we don't die of old age before it's
% finished.
[mfilepath,f,e] = fileparts(mfilename('fullpath'));
cd(mfilepath);

f_profile = false;

imageRoot = flow_imgroot();

% Create new output path for results
% Can assume this doesn't already exist since it's the datetime
outputPath = fullfile(imageRoot, datestr(now, 'flow/yyyymmddTHHMMSS'));
mkdir(outputPath); 

nLevels = 2;
testLevels = 1:nLevels;
nFrames = 60;

tic;
if f_profile, profile clear; profile on; end
[flowPyramidEst, flowConfidence] = ...
    estimate_flow_multilevel(imageRoot, nFrames, outputPath, testLevels);
if f_profile, profile off; profile report; end
toc;

return

% View results
[gt, flowPyramid] = load_ground_truth(imageRoot, nLevels, 1, 0);
allFlows = [];
for level = testLevels(end:-1:1)
    % Save interpolated flowmap.
    flowmap_interp = upscale_flow(flowPyramidEst{level}, level-1);
    filename = sprintf('flowmap_L%02d.png', level);
    imwrite(show_flow_as('rgb', [flowPyramid{1} flowmap_interp]), ...
            fullfile(outputPath, filename));
        
    figure(101); clf;
        b = 8; % border to trim
        imagesc(abs(flowmap_interp(b:end-b,b:end-b)));
        axis('image','off');
        colorbar();
        F = getframe(gcf);
        filename = sprintf('flowmagnitude_L%02d.png', level);
        imwrite(F.cdata, ...
                fullfile(outputPath, filename));
            
    allFlows = [allFlows flowmap_interp nan(size(flowPyramid{1},1),1)];

    % Save interpolated confidence.
    flowconf_interp = interp2(flowConfidence{level}, level-1);
    figure(201); clf;
        b = 8; % border to trim
        imagesc(abs(flowconf_interp(b:end-b,b:end-b)));
        axis('image','off');
        colorbar();
        F = getframe(gcf);
        filename = sprintf('confidence_L%02d.png', level);
        imwrite(F.cdata, ...
                fullfile(outputPath, filename));
end
allFlows = [allFlows flowPyramid{1}];

filename = 'allFlowmaps.png';
imwrite(show_flow_as('rgb', allFlows), ...
        fullfile(outputPath, filename));
    
figure(102); clf;
    imagesc(abs(allFlows));
    axis('image','off');
    colorbar();
    F = getframe(gcf);
    filename = 'allFlowMagnitudes.png';
    imwrite(F.cdata, ...
            fullfile(outputPath, filename));
