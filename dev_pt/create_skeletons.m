clc; clear;

% imtype = 'test';
imtype = 'training';
retpath = [asymmetryroot('shared'),'data\retinograms\DRIVE\',imtype,'\'];

switch imtype
    case 'test', imrng = 1:20;
    case 'training', imrng = 21:40;
end

for i = imrng
    % load masks
    load([retpath,'vessel_masks\',sprintf('%02d_%s_v_mask.mat',i,imtype)]);
    load([retpath,'foveal_masks\',sprintf('%02d_%s_f_mask.mat',i,imtype)]);
    
    % skeletonize vessel mask
    skeleton = bwmorph(vessel_mask, 'thin', 'inf');
    
    % remove vessel and skeleton pixels that are outside the fovea
    vessel_mask(~foveal_mask) = 0;
    skeleton(~foveal_mask) = 0;

    vessel_centres = skeleton(vessel_mask);
    
    % save points list
    filename = sprintf('%02d_%s_centre_idx.mat',i,imtype);
    save([retpath,'vessel_masks\centre_idx_pt\',filename], 'vessel_centres');
end

