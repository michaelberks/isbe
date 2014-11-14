clear all

% profile clear;

forest_job = 'pc20110107T150659';
% image_dir = 'mammograms/2004_screening/abnormals/mat_small/';
image_dir = 'grain_images_small/';

if (1)
% classify a folder full of images
classify_image_set(...
	forest_job, ... % forest job
	image_dir, ... % image dir
	'flip',0);
end

inpath = [asymmetryroot,'data/',image_dir];
d = dir([inpath,'*.mat']);
load([inpath,d(1).name]);

outpath = [asymmetryroot,'data/',image_dir,'results/',forest_job,'/'];
d = dir([outpath,'*.mat']);
out_var = load_uint8([outpath,d(1).name]);

figure(1); clf;
	colormap(gray(256));
	subplot(2,1,1);
		image(mammogram);
		axis('image','ij','off');
	subplot(2,1,2);
		imagesc(out_var);
		axis('image','ij','off');
		
% profile report
