clear;

% convert extended DRIVE images from .mat to .png
imgpath = 'U:\data\retinal\drive\test\images_extended\';
imgdir = dir([imgpath,'*.mat']);

for i = 1:length(imgdir)
	[p,f,e] = fileparts([imgpath,imgdir(i).name]);
	load([p,filesep,f,e]);
	e = '.png';
	imwrite(uint8(ret),[p,filesep,f,e]);
end
