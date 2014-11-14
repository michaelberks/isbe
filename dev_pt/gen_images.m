% parse a directory, load each mat file and save out the mammogram as a
% .png file, scaled to range [0..255]

clear all; clc;

% input path and files
datapath	= 'A:\data\mammograms\2004_screening_processed\mass_roi\';
pathdir		= dir([datapath,'*.mat']);

% output path
pngpath = lower([datapath,'png\']);
pngpath = strrep(pngpath,'a:\','u:\');
if ~exist(pngpath,'dir')
	mkdir(pngpath);
end

for i = 1:length(pathdir);
	str = load([datapath,pathdir(i).name]);
	[p,f,e] = fileparts(pathdir(i).name);
	
	filename = [pngpath,f,'.png'];

	if isfield(str,'mammogram')
% 		mammogram = mammogram-min(mammogram(:));
% 		mammogram = mammogram/( max(mammogram(:))/255 );
		
		imwrite(uint8(str.mammogram),filename);
		
	elseif isfield(str,'mammo_small')
		imwrite(uint8(str.mammo_small),filename);
	
	elseif isfield(str,'mammo_processed')
		imwrite(uint8(str.mammo_processed),filename);

	elseif isfield(str,'out_var')
		imwrite(uint8(str.out_var),filename);
	end
end

