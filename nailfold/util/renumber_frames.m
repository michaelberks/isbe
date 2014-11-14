clc;
clear;
% close all;

datroot = 'U:\projects\nailfold\capture\SRFT\test';
imgpath = fullfile(datroot, '\Left\Digit4\x300\');

outpath = fullfile(imgpath,'original');
if ~exist(outpath,'dir')
    mkdir(outpath);
else
    delete(fullfile(outpath,'*.png'));
end

d = dir(fullfile(imgpath,'*.png'));
for i = 1:length(d)
    outfile = sprintf('frame_%04d.png', i);
    copyfile(fullfile(imgpath,d(i).name), fullfile(outpath,outfile));
end