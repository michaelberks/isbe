function [gt, flowPyramid] = load_ground_truth(imgroot, ...
                                nPyramidLevels, testLevel, ...
                                nFrames)
if (nargin==0 && nargout==0), test(); return; end
                            
gt = struct('flow', [], ...
            'mask', [], ...
            'displacements', [], ...
            'orientation', [], ...
            'presence', []);
        
flowPyramid = cell(1, nPyramidLevels);

gt_filename = fullfile(imgroot, '_ground_truth.mat');
if exist(gt_filename, 'file')
    gt_vals = load(gt_filename, ...
                   'mask', 'jitter', 'flowStackMean');
    
    gt_flow = gt_vals.flowStackMean;
    gt_flow(isnan(gt_flow)) = 0;
    
    flowPyramid = build_flow_pyramid(gt_flow, nPyramidLevels);
    maskPyramid = build_image_pyramid(gt_vals.mask, nPyramidLevels);
    
    gt.flow = flowPyramid{testLevel};
    gt.mask = maskPyramid{testLevel};
    gt.displacements = gt_vals.jitter(2:nFrames,:) / (2^(testLevel-1));
end

pred_path = fullfile(imgroot, 'predictions/orientation');
pred_filename = fullfile(pred_path, 'frame_0001_ori_pred.mat');
if exist(pred_filename, 'file')
    gt_vals = load(pred_filename);
    ori = gt_vals.frame_ori_prediction;
    direction_angle = mod(pi + angle(ori)/2, pi);
    direction = complex(cos(direction_angle), sin(direction_angle));
%     direction = abs(ori) .* direction;
    
    oriPyramid = build_image_pyramid(direction, nPyramidLevels);
    gt.orientation = oriPyramid{testLevel};
end

pred_path = fullfile(imgroot, 'predictions/vessel');
pred_filename = fullfile(pred_path, 'frame_0001_vessel_pred.mat');
if exist(pred_filename, 'file')
    gt_vals = load(pred_filename);
    presPyramid = build_image_pyramid(gt_vals.frame_vessel_prediction, nPyramidLevels);
    gt.presence = presPyramid{testLevel};
end
        

%% Test script
function test()
clc;

imgroot = flow_imgroot();
nPyramidLevels = 4;
testLevel = 2;
nFrames = 60;

gt = load_ground_truth(imgroot, nPyramidLevels, testLevel, nFrames);

figure(5); clf;
    show_flow_as('rgb', gt.orientation);

return
