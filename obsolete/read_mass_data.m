function mass_data = read_mass_data(file_names_a, file_names_o, file_out)

%
% author: Michael Berks
% date: 10/08/2006 10:49
% function: extract the outline and ROI of a mass from an
%           annotated mammogram
%
% inputs:   file_names_a = list of mammogram file names, each mammogram
%                       should contain 1 annotated mass
%           file_names_o = list of mammogram files names that MUST provide
%                       the original mammograms to the annotated list
%           file_out = file path to save the output data to
% outputs:  mass_data = structure storing mass_outline and mass_ROI as
%           returned by GET_SHAPE_MASK for each annotated mass 
%          
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid1 = fopen(file_names_a);
names_a = textscan(fid1, '%s');

fid2 = fopen(file_names_o);
names_o = textscan(fid2, '%s');

clear fid1 fid2 file_names_a file_names_o
fclose('all');

[x, y] = size(names_a{1});
for ii=1:x,
    ii
    %im_anno = imread(names_a{1}{ii});
    %im_orig = imread(names_o{1}{ii});
    [mass_data(ii).shape_border mass_data(ii).shape_ROI]...
        = get_mass_data(names_a{1}{ii}, names_o{1}{ii});
    %clear im_anno im_orig;
end

save(file_out, 'mass_data');