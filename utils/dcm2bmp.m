function dcm2bmp(folder_in, folder_out)
%DCM2BMP converts a folder of DICOM images to bitmap images
%   [] = dcm2bmp(folder_in,folder_in)
%
%   inputs:
%      folder_in - String giving the directory of folder containing DICOM images
%      folder_out- String specifying directory to write bitmaps to
%
%   outputs:
%   
%   notes:
%

%Check folders are filesep terminated
if ~strcmp(folder_in(end), filesep)
    folder_in = [folder_in, filesep];
end
if ~strcmp(folder_out(end), filesep)
    folder_out = [folder_out, filesep];
end

%Get list of 'dcm' files in input directory
dcm_list = dir([folder_in, filesep, '*.dcm']);
k = 0;

%Go through each DCM converting to BMP
for kk = 1:length(dcm_list)
    try
        i1 = dicomread([folder_in, dcm_list(kk).name]);
        i2 = im2uint8(i1); clear i1;
        imwrite(i2, [folder_out dcm_list(kk).name(1:end-3) 'bmp']);
        clear i2;
    catch
        display(['Error converting: ',...
            folder_out dcm_list(kk).name(1:end-3) 'bmp']);
        k = k + 1;
    end            
end

display(['Finished converting ', folder_in, ' with ', num2str(k), ' errors']);