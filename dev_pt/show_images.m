function show_images

clear;

dataroot	= 'a:\data\';
% subgroup	= 'normals';
subgroup	= 'abnormals';
datapath	= [dataroot,'mammograms\2004_screening\',subgroup,'\'];
datapath	= [dataroot,'orientation_maps\rf\2004_screening_processed\',subgroup,'\'];
	ext		= '_class.mat';
% roidir		= dir([datapath,'meta\*CC_meta.mat']);

% scl = 1; sclpath = '';
% scl = 2; sclpath = 'halfsize\';
scl = 4; sclpath = 'qtrsize\';

% datapath	= [datapath,'png\',sclpath];
% 	datadir		= dir([datapath,'*LCC.png']);
datadir		= dir([datapath,'*LCC',ext]);

button = 1;
while (button==1)
 	ip = ceil(rand*length(datadir));
	[p,basename,e] = fileparts(datadir(ip).name);
	basename = basename(1:6);
	display(basename);
	
	if (0)
		% define filenames for left and right sides
		lml_name	= [basename(1:3),'LML.mat'];
		load([datapath,lml_name]); lml = mammogram(1:4:end,1:4:end);
		rml_name	= [basename(1:3),'RML.mat'];
		load([datapath,rml_name]); rml = mammogram(1:4:end,1:4:end);
	else
		% define filenames for left and right sides
		lml_name	= [datapath,basename(1:3),'LCC',ext];
% 		lml = imread(lml_name); lml = lml(:,:,1);
		lml = load_uint8(lml_name,true);
			lml_map = abs(lml(1:scl:end,end:-scl:1,1)); % flipped l-r
			lml = angle(conj(lml(1:scl:end,end:-scl:1,1))); % flip angle, too
			lml = (lml+pi/2)/pi;

		rml_name	= [datapath,basename(1:3),'RCC',ext];
% 		rml = imread(rml_name); rml = rml(:,:,1);
		rml = load_uint8(rml_name,true);
			rml_map = abs(rml(1:scl:end,1:scl:end,1));
			rml = angle(rml(1:scl:end,1:scl:end,1));
			rml = (rml+pi/2)/pi;
	end

	load([dataroot,'masks\2004_screening\',subgroup,'\',basename(1:3),'LCC_mask.mat']); 
		maskl = mask(1:scl:end,end:-scl:1); lml(~maskl) = NaN;
	load([dataroot,'masks\2004_screening\',subgroup,'\',basename(1:3),'RCC_mask.mat']);
		maskr = mask(1:scl:end,1:scl:end); rml(~maskr) = NaN;

% 	lml = 128+64*normim(lml,'localstd',11,5);
% 	rml = 128+64*normim(rml,'localstd',11,5);

	trim = 20;
	lml = lml(trim+1:end-trim,trim+1:end-trim);
	lml_map = lml_map(trim+1:end-trim,trim+1:end-trim);
	rml = rml(trim+1:end-trim,trim+1:end-trim);
	rml_map = rml_map(trim+1:end-trim,trim+1:end-trim);

% 	% flip left and right views to avoid familiarity
% 	f_fliplr = double(rand>0.5);
% 	if (f_fliplr)
% 		s1 = 2; s2 = 1;
% 		lml = lml(:,end:-1:1);
% 		rml = rml(:,end:-1:1);
% 	else
		s1 = 1; s2 = 2;
% 	end

% 	% invert image intensities?
% 	f_flipint = (rand>0.5);
% 	if (f_flipint)
% 		rml = 255-rml; lml = 255-lml;
% 	end
	
	lml_rgb = lml; lml_rgb(:,:,2) = 1; lml_rgb(:,:,3) = lml_map; lml_rgb = hsv2rgb(lml_rgb);
	rml_rgb = rml; rml_rgb(:,:,2) = 1; rml_rgb(:,:,3) = rml_map; rml_rgb = hsv2rgb(rml_rgb);

	figure(1); clf;
		subplot(2,2,s1); image(lml_rgb); title(''); axis('off','image');
		subplot(2,2,s2); image(rml_rgb); title(''); axis('off','image');
% 	figure(1); clf; colormap([0 0 0; hsv(255)]);
% 		subplot(1,2,s1); image(uint8(255*lml)); title(''); axis('off','image');
% 		subplot(1,2,s2); image(uint8(255*rml)); title(''); axis('off','image');
	
	subplot(2,2,s1); [x,y] = ginput(1);
	
	r_min = 10; r_max = 80;
	rad = round(logspace(log10(r_min),log10(r_max),5))



	pause; continue;
	
	
	
	load([datapath,'meta\',roidir(ip).name]);
	meta_xy = meta_xy/scl;
	if (f_fliplr), meta_xy(:,1) = size(lml,2)-meta_xy(:,1); end
	
	[x,y,button] = ginput(1);
	if (basename(4)=='R')
		subplot(1,2,s1); title('abnormal');
			hold on; plot(meta_xy([1:end,1],1),meta_xy([1:end,1],2),'r-');
		subplot(1,2,s2); title('normal');
		[x,y,button] = ginput(1);
		subplot(1,2,s1); hold off; imshow(rml); title('abnormal');
	else
		subplot(1,2,s1); title('normal');
		subplot(1,2,s2); title('abnormal');
			hold on; plot(meta_xy([1:end,1],1),meta_xy([1:end,1],2),'r-');
		[x,y,button] = ginput(1);
		subplot(1,2,s2); hold off; imshow(lml); title('abnormal');
	end		
	[x,y,button] = ginput(1);
end


function [lml,rml,lml_map,rml_map] = ...
	load_images(dataroot,datapath,subgroup,basename,ext,scl)

