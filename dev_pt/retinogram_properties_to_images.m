clc; clear;

imtype = 'test';
% imtype = 'training';
retpath = [asymmetryroot('shared'),'data\retinograms\DRIVE\',imtype,'\'];

load([retpath,'retinogram_properties.mat'],...
	 'vessel_inds','centre_inds','line_widths','line_contrasts');

switch imtype
    case 'test', imrng = 1:20;
    case 'training', imrng = 21:40;
end

% create folders where necessary
width_path = [retpath,'width_maps/'];
if ~exist(width_path,'dir')
    mkdir(width_path);
end
contrast_path = [retpath,'contrast_maps/'];
if ~exist(contrast_path,'dir')
    mkdir(contrast_path);
end

for i = imrng
    % load masks
    load([retpath,'vessel_masks\',sprintf('%02d_%s_v_mask.mat',i,imtype)]);
    load([retpath,'foveal_masks\',sprintf('%02d_%s_f_mask.mat',i,imtype)]);
    
    % mask out vessel points that are outside the foveal region
    vessel_mask(~foveal_mask) = 0;
    
    % deal with NaNs
    load([retpath,'orientations\',sprintf('%02d_ori1.mat',i)]); % gt_ori
    gt_ori(~vessel_mask) = NaN;
	vessel_mask = vessel_mask & ~isnan(gt_ori);

    % line widths are maps of the approximate width of the vessel at a given
    % point.
    width_map = double(vessel_mask);
    width_map(width_map==1) = line_widths{i};
    save([width_path,sprintf('%02d_%s_width.mat',i,imtype)], 'width_map');
    imwrite(uint8(width_map/max(width_map(:)) * 255), ...
            [width_path,sprintf('%02d_%s_width.png',i,imtype)]);

    % Line contrasts gives the approximate standard deviation within a 15x15
    % window centred at the pixel. (The window size is debatable.)
    contrast_map = double(vessel_mask);
    contrast_map(contrast_map==1) = line_contrasts{i};
    save([contrast_path,sprintf('%02d_%s_contrast.mat',i,imtype)], 'contrast_map');
    imwrite(uint8(contrast_map/max(contrast_map(:)) * 255), ...
            [contrast_path,sprintf('%02d_%s_contrast.png',i,imtype)]);
    
end
 