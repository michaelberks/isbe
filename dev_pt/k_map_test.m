clc; clear;

% test = '002LCC'; mass_c2 = [320,378]; 
test = '013LML'; mass_c2 = [400,400]; 
% test = '024RCC'; mass_c2 = [400,400]; 
% test = '024RML'; mass_c2 = [377,393]; 
% test = '029RCC'; mass_c2 = [400,400]; 
% test = '039RML'; mass_c2 = [460,400]; 
% test = '046RCC'; mass_c2 = [385,340]; 
% test = '056RML'; mass_c2 = [400,400]; 
% test = '060RML'; mass_c2 = [355,425]; mass_c2 = [560,245];
% test = '061RCC'; mass_c2 = [520,535]; % false postive?
% test = '073LML'; mass_c2 = [400,400]; 
% test = '085LCC'; mass_c2 = [400,400]; 
% test = '085LML'; mass_c2 = [400,400]; 
% test = '097RML'; mass_c2 = [400,380]; 
% test = '101RCC'; mass_c2 = [400,400]; 
% test = '112LCC'; mass_c2 = [400,440]; 

figure(1); clf; colormap(gray(256));
inpath	= ['A:\data\mammograms\2004_screening_processed\mass_roi\'];
img		= load_uint8([inpath,test,'_roi.mat']);
	imagesc(img); axis('image');
inpath	= ['A:\data\mammograms\2004_screening\abnormals\'];
load([inpath,'meta\',test,'_meta.mat']);
mass_c	= round(mean(meta_xy));


cmap	= [0 0 0; hsv(255)];
figure(2); clf; colormap(cmap); sbsz = [2,3];

legbar	= ones(20,1)*linspace(0,255,800);

methods	= {'g2d','rf'};
for i_method = 1:length(methods)
	method	= methods{i_method};
	
	% orientation map
	inpath	= ['A:\data\orientation_maps\',method,'\2004_screening_processed\mass_roi\'];

	% line map
	switch method
		case 'g2d',
			orimap_g2d	= load_uint8([inpath,test,'_roi.mat']);
			
			inpath	= ['A:\data\line_maps\',method,'\2004_screening_processed\mass_roi\'];
			load([inpath,test,'_roi.mat']); 
			linemap_g2d = in_var;
		case 'rf',
			orimap_rf = load_uint8([inpath,test,'_roi.mat']);
			linemap_rf = abs(orimap_rf);
			orimap_rf = mod(angle(orimap_rf),pi);
			
			linemap_rf = (linemap_rf>0.5);
	end
end

% choose linemap
linemap = linemap_g2d;
% linemap = linemap_rf;

for i_method = 1:length(methods)
	method	= methods{i_method};
	
	switch method
		case 'g2d',	orimap = orimap_g2d;
		case 'rf',	orimap = orimap_rf;
	end
	
	if 0
		% load computed k-map
		inpath	= ['A:\data\k_stellate_maps\',method,'\2004_screening_processed\abnormals\'];
		load([inpath,test,'_mask.mat']);
		k_data	= load_uint8([inpath,test,'_f1.mat']);
		k_map	= nan(size(mask));
		k_map(mask) = k_data(:,1);
		k_map	= k_map(mass_c(2)-399:mass_c(2)+400,mass_c(1)-399:mass_c(1)+400);
	else
		map_sz = 64; spacing = 4;
		r_max = 200; r_min = r_max/2; R = r_min/4;
		[f1,f2,k_mask] = karssemeijer_radial_projection_multiscale(...
							linemap,orimap,...
							'num_angles',4,...
						    'r_min',r_min,'r_max',r_max,'R',R,...
							'x_min',mass_c2(1)-map_sz,...
							'x_max',mass_c2(1)+map_sz,...
							'y_min',mass_c2(2)-map_sz,...
							'y_max',mass_c2(2)+map_sz,...
							'spacing',spacing);
% 		[f1 f2]
		f1		= 0.5 + min(max(f1,-100),100)/200;
		k_map	= nan(size(k_mask));
		k_map(k_mask) = f1(:,1);
		k_map	= k_map(mass_c2(2)-map_sz:mass_c2(2)+map_sz,...
						mass_c2(1)-map_sz:mass_c2(1)+map_sz);
					
		width = sqrt(length(f1));
		k_map = reshape(f1,[width,width]);
	end
	
	% convert orimap to hsv format
	orimap_hsv(:,:,1) = [orimap/pi; ones(20,1)*linspace(0,1,size(orimap,2))]; % hue
	orimap_hsv(:,:,2) = 1; % saturation
	orimap_hsv(:,:,3) = [linemap; ones(20,size(orimap,2))]; % value
	orimap_rgb = hsv2rgb(orimap_hsv);
	
	figure(1); hold on;
	hold on; rectangle('position',[mass_c2-r_max,r_max*2,r_max*2],'edgecolor','g','curvature',[1 1]);
	hold on; rectangle('position',[mass_c2-r_min,r_min*2,r_min*2],'edgecolor','g','curvature',[1 1]);
	hold on; rectangle('position',[mass_c2-R,R*2,R*2],'edgecolor','b','curvature',[1 1]);
	if map_sz>0
		hold on; rectangle('position',[mass_c2-map_sz,map_sz*2,map_sz*2],'edgecolor','r');
	else
		hold on; plot(mass_c2(1),mass_c2(2),'r+');
	end
		
	figure(2);
	mysubplot(sbsz); image(uint8(255*linemap(:,:,[1,1,1]))); axis('image','off');
	if map_sz>0
		hold on; rectangle('position',[mass_c2-map_sz,map_sz*2,map_sz*2],'edgecolor','r');
	else
		hold on; plot(mass_c2(1),mass_c2(2),'r+');
	end
% 	mysubplot(sbsz); image(uint8(1+255*orimap/pi)); axis('image','off');
	mysubplot(sbsz); image(orimap_rgb); axis('image','off');
	mysubplot(sbsz); image(uint8(1+255*k_map)); axis('image','off'); colormap(jet(256));
	drawnow;
end

