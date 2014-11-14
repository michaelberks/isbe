clear;

froot = 'U:\projects\mammography\data\retinograms\drive\test\predictions\';
% fpath = [froot,'centre_detection\dt\rf_3\'];
% fpath = [froot,'detection\dt\rf_3\'];
% fpath = [froot,'junction_centre_detection\dt\rf_1\'];
fpath = [froot,'junction_detection\dt\rf_3\'];

if ~exist([fpath,'png'],'dir')
	mkdir([fpath,'png']);
end

fdir = dir([fpath,'*.mat']);
for i = 1:length(fdir)
	img = load_uint8([fpath,fdir(i).name]);
	if (all(0 <= img(:)) && all(img(:) <= 1))
		% probability map - do nothing
	else
		img = img - min(img(:));
		img = img / max(img(:));
	end
	imwrite(uint8(img*255), [fpath,'png\',strrep(fdir(i).name,'.mat','.png')]);
end