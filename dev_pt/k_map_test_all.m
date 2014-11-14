% clc;
clear;

% names of methods
method_list = {'g2d','rf_thin'};

% mask
maskpath = ['A:\data\masks\2004_screening\abnormals\'];
maskdir = dir([maskpath,'*_mask.mat']);

for iimg = 1:length(maskdir)
	% get basename (e.g. 045RML)
	basename = maskdir(iimg).name(1:6)
	
	% load mask for this image
	load([maskpath,basename,'_mask.mat']);

	% compute map for both methods
	for imethod = 1:length(method_list)
		method	= method_list{imethod};

		% paths
		oripath	= ['A:\data\orientation_maps\',method,'\2004_screening_processed\abnormals\'];
		linepath = ['A:\data\line_maps\',method,'\2004_screening_processed\abnormals\'];
		outpath = ['U:\projects\mammography\data\k_stellate_maps\'];

		% skip files that have been processed
		if exist([outpath,basename,'_',method,'_f1.png'],'file')
			continue;
		end
		
		% line map
		switch method
			case 'g2d',
				orimap	= load_uint8([oripath,basename,'_ori.mat']);
				load([linepath,basename,'_lines.mat']); 
				linemap = line_map & mask;
				clear line_map;
				
			case {'rf','rf_thin'}
				orimap = load_uint8([oripath,basename,'_class.mat']);
				linemap = abs(orimap);
		% 		linemap = load_uint8([inpath,test,'_class.mat']); 
				orimap = mod(angle(orimap),pi);
				linemap = (linemap>0.5) & mask;
				linemap = linemap .* mask;
		end

		% compute k-map
		spacing = 16; map_sz = 8*spacing;
		r_max = 180; r_min = r_max/2; R = r_min/2;
		% r_max = 180; r_min = 10; R = 5;

		x_min = 1; x_max = size(orimap,2); 
		y_min = 1; y_max = size(orimap,1);

		tic;
		[f1,f2,k_mask] = karssemeijer_radial_projection_multiscale(...
							linemap,orimap,...
							'num_angles',24,...
							'x_min',x_min,'x_max',x_max,...
							'y_min',y_min,'y_max',y_max,...
							'r_min',r_min,'r_max',r_max,'R',R,...
							'spacing',spacing,...
							'mask',mask);
		toc;

		% copy values to images and subsample used points
		k_map_f1 = nan(size(k_mask));
			k_map_f1(k_mask) = f1(:,1);
			k_map_f1 = k_map_f1(y_min:spacing:y_max,x_min:spacing:x_max);
		k_map_f2 = nan(size(k_mask));
			k_map_f2(k_mask) = f2(:,1);
			k_map_f2 = k_map_f2(y_min:spacing:y_max,x_min:spacing:x_max);

		% interpolate to get a decent sized image on the cheap
		tic;
		k_map_f1_interp = interp2(k_map_f1,2);
		k_map_f2_interp = interp2(k_map_f2,2);
		toc;
			
		% save a copy of the positive parts of f1
		figure(4); colormap(jet(256));
			k_map_f1_pos = k_map_f1_interp;
			k_map_f1_pos(k_map_f1_pos<=0) = nan;
			imagesc(k_map_f1_pos); axis('off','image');

		save([outpath,basename,'_',method,'_kmaps.mat'],...
			'k_map_f1','k_map_f2');
		imwrite(uint8(255*normim(k_map_f1_pos)),jet(256),...
			[outpath,basename,'_',method,'_f1.png']);

	end % for imethod
end % for iimg
