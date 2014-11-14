% clc;
clear;

% test = '002LCC'; tp = [774,1246]; %fp = [830,1540]; fp = [542,1214]; fp = [302,996];
% test = '013LML'; %fp = [636,1928];
% test = '024RCC'; 
% test = '024RML'; 
% test = '029RCC'; 
% test = '039RML'; 
test = '045LML'; tp = [508,1672];
% test = '046RCC'; 
% test = '056RML'; 
% test = '060RML'; 
% test = '061RCC';
% test = '073LML';
% test = '085LCC';
% test = '085LML';
% test = '097RML';
% test = '101RCC';
% test = '112LCC';

method	= 'g2d';
% method	= 'rf_thin';

% mask
inpath	= ['A:\data\masks\2004_screening\abnormals\'];
load([inpath,test,'_mask.mat']);

% get metadata
inpath	= ['A:\data\mammograms\2004_screening\abnormals\'];
load([inpath,'meta\',test,'_meta.mat']);
clear meta_xy;
if ~exist('tp','var'); tp = round(mean(meta_xy)); end % true positive

% orientation map
inpath	= ['A:\data\orientation_maps\',method,'\2004_screening_processed\abnormals\'];

% line map
switch method
	case 'g2d',
		orimap	= load_uint8([inpath,test,'_ori.mat']);
		inpath	= ['A:\data\line_maps\',method,'\2004_screening_processed\abnormals\'];
		load([inpath,test,'_lines.mat']); 
		linemap = line_map & mask;
		clear line_map;
	case {'rf','rf_thin'}
		orimap = load_uint8([inpath,test,'_class.mat']);
		linemap = abs(orimap);
% 		inpath	= ['A:\data\line_maps\old_rf\2004_screening_processed\abnormals\'];
% 		linemap = load_uint8([inpath,test,'_class.mat']); 
		orimap = mod(angle(orimap),pi);
		linemap = (linemap>0.5) & mask;
		linemap = linemap .* mask;
end

% compute k-map
spacing = 32; map_sz = 8*spacing;
r_max = 180; r_min = r_max/2; R = r_min/2;
% r_max = 180; r_min = 10; R = 5;

p = tp; % test true positive
if exist('fp','var'), p = fp; end % test false positive

x_min = max(1,p(1)-map_sz); x_max = min(p(1)+map_sz,size(orimap,2));
y_min = max(1,p(2)-map_sz); y_max = min(p(2)+map_sz,size(orimap,1));

x_min = 1; x_max = size(orimap,2); y_min = 1; y_max = size(orimap,1);

tic;
[f1vec,f2vec,k_mask] = karssemeijer_radial_projection_multiscale(...
					linemap,orimap,...
					'num_angles',1,...
					'x_min',x_min,'x_max',x_max,...
					'y_min',y_min,'y_max',y_max,...
					'r_min',r_min,'r_max',r_max,'R',R,...
					'spacing',spacing,...
					'mask',mask);
toc;
f1	= 0.5 + min(max(f1vec,-50),50)/100;
f2	= 0.5 + min(max(f2vec,-4),4)/8;

display([min(f1vec) max(f1vec)]);

k_map_f1 = nan(size(k_mask));
	k_map_f1(k_mask) = f1(:,1);
	k_map_f1 = k_map_f1(y_min:spacing:y_max,x_min:spacing:x_max);
% 	k_map_f1 = normim(k_map_f1);
k_map_f2 = nan(size(k_mask));
	k_map_f2(k_mask) = f2(:,1);
	k_map_f2 = k_map_f2(y_min:spacing:y_max,x_min:spacing:x_max);
% 	k_map_f2 = normim(k_map_f2);

% interpolate to get a decent sized image on the cheap
tic;
	k_map_f1_interp = interp2(k_map_f1,1);
	k_map_f2_interp = interp2(k_map_f2,1);
toc;
		
% inpath	= ['A:\data\mammograms\2004_screening_processed\abnormals\'];
% load([inpath,test,'.mat']); img = mammo_processed; clear mammo_processed;
inpath	= ['A:\data\mammograms\2004_screening\abnormals\'];
load([inpath,test,'.mat']); img = mammogram; clear mammogram;

% convert orimap to hsv format
clear orimap_3;
orimap_3(:,:,1) = orimap(1:2:end,1:2:end)/pi; % hue
orimap_3(:,:,2) = 1; % saturation
orimap_3(:,:,3) = linemap(1:2:end,1:2:end); % value
orimap_3 = hsv2rgb(orimap_3);

figure(1); clf; spsz = [2,3]; colormap(redgreen(256));
	mysubplot(spsz,1); image(uint8(repmat(img,[1,1,3]))); axis('off','image'); hold on; 
		plot(p(1),p(2),'b.');
		if map_sz>0
			rectangle('position',[p-map_sz,2*map_sz,2*map_sz],'edgecolor','r');
		end
	mysubplot(spsz,2); imagesc(repmat(linemap,[1,1,3])); axis('off','image');
	mysubplot(spsz,3); image(orimap_3); axis('off','image');
	
	[xx,yy] = meshgrid(-r_max:r_max,-r_max:r_max); rr = sqrt(xx.*xx+yy.*yy);
	subimg = double(img(p(2)-r_max:p(2)+r_max,p(1)-r_max:p(1)+r_max));
% 	subimg = normim(subimg,'localstd',3,40);
	subimg((rr>r_max) | (rr<r_min)) = NaN;
	mysubplot(spsz,4); image(uint8(255*repmat(normim(subimg),[1,1,3]))); axis('off','image'); hold on;
% 		rectangle('position',[0*[1,1],2*r_max,2*r_max],'edgecolor','r','curvature',[1 1]);
% 		rectangle('position',[(r_max-r_min)*[1,1],2*r_min,2*r_min],'edgecolor','r','curvature',[1 1]);
% 		rectangle('position',[(r_max-R)*[1,1],2*R,2*R],'edgecolor','r','curvature',[1 1]);
	mysubplot(spsz,5); image(uint8(1+255*k_map_f1_interp)); axis('off','image');
		hold on; plot(meta_xy([1:end,1],1)/spacing,meta_xy([1:end,1],2)/spacing,'k-');
		colorbar('East');
	mysubplot(spsz,6); image(uint8(1+255*k_map_f2_interp)); axis('off','image');
		hold on; plot(meta_xy([1:end,1],1),meta_xy([1:end,1],2),'k-');
		colorbar('East');

figure(3);
	mysubplot([3,2],5); image(uint8(255*repmat(normim(subimg),[1,1,3]))); axis('off','image'); hold on;
% 		rectangle('position',[(r_max-R)*[1,1],2*R,2*R],'edgecolor','r','curvature',[1 1]);

figure(4); colormap(jet(256));
	k_map_f1_pos = k_map_f1_interp-0.5;
% 	k_map_f1_pos(k_map_f1_pos<=0) = nan;
% 	k_map_f1_pos(isnan(k_map_f1_pos)) = 0;
	imagesc(k_map_f1_pos); axis('off','image');
	if exist('meta_xy','var')
		hold on; plot(meta_xy([1:end,1],1)/spacing,...
					  meta_xy([1:end,1],2)/spacing,'k-');
	end



orimap_roi = double(orimap(p(2)-r_max:p(2)+r_max,p(1)-r_max:p(1)+r_max));
	orimap_roi((rr>r_max) | (rr<r_min)) = NaN;
linemap_roi = double(linemap(p(2)-r_max:p(2)+r_max,p(1)-r_max:p(1)+r_max));
	linemap_roi((rr>r_max) | (rr<r_min)) = NaN;

clear tmpimg;
tmpimg(:,:,1) = orimap_roi/pi;
tmpimg(:,:,2) = 1;
tmpimg(:,:,3) = linemap_roi;

% figure(10); image(uint8(255*repmat(normim(subimg),[1,1,3]))); axis('off','image');
% figure(12); image(uint8(255*linemap_roi)); axis('off','image'); colormap([0 0 0; redgreen(255)]);
% figure(11); image(hsv2rgb(tmpimg)); axis('off','image');

s = sort(f1vec);
v = linspace(0,1,length(s));
i = find(diff(s>0)); % zero crossing
figure(10); plot(s,v,'b-',[0,0],[0,1],'k:',[s(1),s(end)],v([i,i]),'k:');

