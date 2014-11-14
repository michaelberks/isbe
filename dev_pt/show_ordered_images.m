clear all

datapath	= 'a:\data\mammograms\2004_screening\abnormals\';
roidir		= dir([datapath,'meta\*ML_meta.mat']);

% scl = 1; sclpath = '';
% scl = 2; sclpath = 'halfsize\';
scl = 4; sclpath = 'qtrsize\';

r_mouse = 3;

button = 1;
for ip = 1:length(roidir)
	[p,basename,e] = fileparts(roidir(ip).name);
	set(gcf,'name',basename);

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
% 		lcc_name	= [datapath,'png\',sclpath,basename(1:3),'LCC.png'];
% 		lcc = imread(lcc_name); lcc = lcc(:,:,1)';
% 		rcc_name	= [datapath,'png\',sclpath,basename(1:3),'RCC.png'];
% 		rcc = imread(rcc_name); rcc = rcc(:,:,1)';
	end

	sp_size = [1,2];

	figure(1); clf;
		mysubplot(sp_size,1); imshow(rml); title('');
		mysubplot(sp_size,2); imshow(lml); title('');
% 		mysubplot(sp_size,3); imshow(rcc); title('');
% 		mysubplot(sp_size,4); imshow(lcc); title('');

	load([datapath,'meta\',roidir(ip).name]);	
	meta_ml = meta_xy(:,[1,2])/scl;
% 	load([datapath,'meta\',strrep(roidir(ip).name,'ML','CC')]);
% 	meta_cc = meta_xy(:,[2,1])/scl;
		
	if (basename(4)=='R')
		mysubplot(sp_size,2); title('normal');
		mysubplot(sp_size,1); title('abnormal');
		[x,y,button] = ginput(1); if (button == 3), break; end;
			mysubplot(sp_size,1); hold on; plot(meta_ml([1:end,1],1),meta_ml([1:end,1],2),'r-');
% 			mysubplot(sp_size,3); hold on; plot(meta_cc([1:end,1],1),meta_cc([1:end,1],2),'r-');
		[x,y,button] = ginput(1); if (button == 3), break; end;
		mysubplot(sp_size,1); hold off; imshow(rml); title('abnormal');
% 		mysubplot(sp_size,3); hold off; imshow(rcc);
	else
		mysubplot(sp_size,1); title('normal');
		mysubplot(sp_size,2); title('abnormal');
		[x,y,button] = ginput(1); if (button == 3), break; end;
			mysubplot(sp_size,2); hold on; plot(meta_ml([1:end,1],1),meta_ml([1:end,1],2),'r-');
% 			mysubplot(sp_size,4); hold on; plot(meta_cc([1:end,1],1),meta_cc([1:end,1],2),'r-');
		[x,y,button] = ginput(1); if (button == 3), break; end;
		mysubplot(sp_size,2); hold off; imshow(lml); title('abnormal');
% 		mysubplot(sp_size,4); hold off; imshow(lcc);
	end		
	[x,y,button] = ginput(1); if (button == 3), break; end;
end

close(gcf);