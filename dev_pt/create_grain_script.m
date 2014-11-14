clc;

datapath	= [asymmetryroot,'data\grain_images\'];
pngpath		= [datapath,'png\'];
if ~exist(pngpath,'dir'), mkdir(pngpath); end

for i = 1:100
	[mammogram,thetas(i)] = gen_grain_image;
	
	% define basename for files
	filebase = sprintf('%03d',i);

	% save image to .mat file in datapath
	save([datapath,filebase,'.mat'],'mammogram');
	
	% save image to .png file in pngpath
	imwrite(uint8(mammogram),[pngpath,filebase,'.png']);
end

% save rotation angles to text file
fid = fopen([datapath,'rotations.txt'],'w');
	fprintf(fid,'%.2f\n',theta);
fclose(fid);