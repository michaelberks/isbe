clear all

datapath	= 'a:\data\mammograms\2004_screening\abnormals\';
roidir		= dir([datapath,'meta\*ML_meta.mat']);

% scl = 1; sclpath = '';
scl = 2; sclpath = 'halfsize\';
% scl = 4; sclpath = 'qtrsize\';

button = 1;
while (button==1)
	ip = ceil(rand*length(roidir));
	[p,basename,e] = fileparts(roidir(ip).name);

	if (0)
		% define filenames for left and right sides
		lml_name	= [basename(1:3),'LML.mat'];
		load([datapath,lml_name]); lml = mammogram(1:4:end,1:4:end);
		rml_name	= [basename(1:3),'RML.mat'];
		load([datapath,rml_name]); rml = mammogram(1:4:end,1:4:end);
	else
		% define filenames for left and right sides
		lml_name	= [datapath,'png\',sclpath,basename(1:3),'LML.png'];
		lml = imread(lml_name); lml = lml(:,:,1);
		rml_name	= [datapath,'png\',sclpath,basename(1:3),'RML.png'];
		rml = imread(rml_name); rml = rml(:,:,1);
	end

	% flip left and right views to avoid familiarity
	f_fliplr = double(rand>0.5);
	if (f_fliplr)
		s1 = 2; s2 = 1;
		lml = lml(:,end:-1:1);
		rml = rml(:,end:-1:1);
	else
		s1 = 1; s2 = 2;
	end

	% invert image intensities?
	f_flipint = (rand>0.5);
	if (f_flipint)
		rml = 255-rml; lml = 255-lml;
	end
	
	figure(1); clf;
		subplot(1,2,s1); imshow(rml); title('');
		subplot(1,2,s2); imshow(lml); title('');

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

