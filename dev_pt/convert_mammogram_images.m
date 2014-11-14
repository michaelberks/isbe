clear;

% convert extended DRIVE images from .mat to .png
imgpath = 'synthetic_lines\real512\masks\';
inroot	= [asymmetryroot,'data\'];
inpath	= [inroot,imgpath];
outroot	= [asymmetryroot,'data\'];
outpath	= [outroot,imgpath];
imgdir	= dir([inpath,'*.mat']);

if ~exist(outpath,'dir')
	mkdir(outpath);
end

for i = 1:length(imgdir)
	[p,f,e] = fileparts([inpath,imgdir(i).name]);
	infile = fullfile(inpath,[f,'.mat']);
	outfile = fullfile(outpath,[f,'.png']);
	vars = whos('-file',infile);
	var_names = {vars(:).name};
	
	if any(strcmp(var_names,'test_image'))
		load(infile,'test_image');
		imwrite(uint8(test_image),outfile);
	elseif any(strcmp(var_names,'bar_real'))
		load(infile,'bar_real');
		imwrite(uint8(bar_real),outfile);
	elseif any(strcmp(var_names,'mask'))
		load(infile,'mask');
		imwrite(uint8(255*mask),outfile);
	end
end
