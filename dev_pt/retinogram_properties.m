clc; clear;

% imtype = 'test';
imtype = 'training';
retpath = [asymmetryroot('shared'),'data\retinograms\DRIVE\',imtype,'\'];

vessel_inds = cell(1,20);
centre_inds = cell(1,20);
line_widths = cell(1,20);
line_contrasts = cell(1,20);

switch imtype
    case 'test', imrng = 1:20;
    case 'training', imrng = 21:40;
end

for i = imrng
    % load masks
    load([retpath,'vessel_masks\',sprintf('%02d_%s_v_mask.mat',i,imtype)]);
    load([retpath,'foveal_masks\',sprintf('%02d_%s_f_mask.mat',i,imtype)]);
    load([retpath,'vessel_masks\centre_idx\',sprintf('%02d_%s_centre_idx.mat',i,imtype)]);
    
    % sort out masks
    vessel_mask(~foveal_mask) = 0;
    centre_mask = vessel_mask;
    centre_mask(centre_mask==1) = vessel_centres;

	% deal with NaNs
    load([retpath,'orientations\',sprintf('%02d_ori1.mat',i)]); % gt_ori
    gt_ori(~vessel_mask) = NaN;
	vessel_mask = vessel_mask & ~isnan(gt_ori);
    centre_mask = centre_mask & ~isnan(gt_ori);
	
	vessel_inds{i} = find(vessel_mask);
	centre_inds{i} = (centre_mask(vessel_mask));
	
	
    % compute line widths
	img_dt = distance_transform(1e6*vessel_mask);
	if 0
		% old method: good for centrelines but rubbish off-centre
		line_widths{i} = img_dt(vessel_mask);
	else
		% new method: replace every pixel with twice the width assigned to the
		% nearest centreline pixel
		width_img = zeros(size(img_dt));
		[yv,xv] = find(vessel_mask);
		[yc,xc] = find(centre_mask);
		for iv = 1:length(yv)
			dx = xc-xv(iv); dy = yc-yv(iv);		
			[ignore,ind] = min(dx.*dx + dy.*dy);
			width_img(yv(iv),xv(iv)) = 2*img_dt(yc(ind),xc(ind))-1;
		end
		line_widths{i} = width_img(vessel_mask);
	end
	
	
	% compute line contrasts
% 	im = imread([retpath,'images_extended\',sprintf('%02d_test.tif',i)]);
% 	im = im(:,:,2);
	load([retpath,'images_extended\',sprintf('%02d_%s_ext.mat',i,imtype)]);
	im = rgb2gray(ret);
	
	if 0
		% the old way: approximate contrast with local standard deviation
		[mn,sd] = local_image_stats(im,7);
		line_contrasts{i} = sd(vessel_mask);
	else
		% the new way: smooth the image using only background pixels and
		% compute difference
		bg_mask = (foveal_mask & ~vessel_mask);
		[mn,sd] = local_image_stats(im,15,0,bg_mask);
		diff_im = abs(double(im)-mn);
		line_contrasts{i} = diff_im(vessel_mask);
	end	
end

save([retpath,'retinogram_properties.mat'],...
	 'vessel_inds','centre_inds','line_widths','line_contrasts');
 