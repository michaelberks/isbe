function convert_dicom_filenames(dir_in, dir_out)
%CONVERT_DICOM_FILENAMES Convert DICOM filenames from generic names to
%   [file_out] = convert_dicom_filenames(file_in)
%
% Inputs:
%      file_in - Filename of dicom image to change (assumed to be in
%      generic form)
%
%
% Outputs:
%      file_out - Filename modified to inlcude image specific information
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 06-Oct-2010
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester
if nargin < 2
    dir_out = dir_in;
end

%Make sure dir and dir_out are filesep ended
if ~strcmpi(dir_in(end), '\')
    dir_in = [dir_in '\'];
end
if ~strcmpi(dir_out(end), '\')
    dir_out = [dir_out '\'];
end

%Get path of all sub-directories
dir_path = genpath(dir_in);

%Workout where each sub-directory starts and ends
dir_ends = strfind(dir_path, pathsep);
num_dirs = length(dir_ends);

%Loop through the sub-directories looking for dicom files
id2 = -1;
for ii = 1:num_dirs
    id1 = id2 + 2;
    id2 = dir_ends(ii) - 1;
    dir_name = [dir_path(id1:id2) filesep];
    dicom_list = dir([dir_name '*.dcm']);
    
    %For each dicom file found, try converting it's name and moving to the
    %output directory
    for jj = 1:length(dicom_list)
        try
            file_out = convert_dicom_filename([dir_name dicom_list(jj).name], dir_out);
            display(['Successfully renamed ' dir_name dicom_list(jj).name ' to ' dir_out file_out]);
        catch
            display(['Problem renaming ' dir_name dicom_list(jj).name]);
            err = lasterror;
            display(err.message);
        end
    end
end

function [file_out] = convert_dicom_filename(file_in, dir_out)
            
%read in info from dicom header
file_info = dicominfo(file_in);

%display info (comment this out later when we know what fields we need)
%display(file_info);

% %Get relevant data from each field, for example:
pma = file_info.AccessionNumber;
leftright = file_info.ImageLaterality;
view = file_info.ViewPosition;

if ~isempty(strfind(file_info.PresentationIntentType, 'PROCESSING'))
    %Image is FOR PROCESSING i.e. raw
    type = 'RAW';
else
    %Image must be 'FOR PRESENTATION', i.e. processed
    type = 'PRO';
end

% %Combine info into a new file name as desired
file_out = [pma '_' leftright view '_' type '.dcm'];
% 
% %Move the old filename to the new filename
movefile(file_in, [dir_out file_out]);

