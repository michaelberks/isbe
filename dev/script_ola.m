%**************************************************************************
% Script for Ola

mam_dir = 'C:\isbe\mammograms\new_CAD\bmP_2004\';
mam_list = dir('C:\isbe\mammograms\new_CAD\bmP_2004\*.bmp');
%

for ii = 1:length(mam_list);
    mlo = ~isempty(strfind(mam_list(ii).name, 'ML'));
    
    %read image in
    mam = imread([mam_dir, mam_list(ii).name]);
    [height width] = size(mam);
    
    %resize to 1024 rows
    mam = imresize(mam, [1024 NaN], 'bilinear');
    %convert from uint8 to double
    mam = double(mam);
    %if it a right mam flip so chest wall on left
    if ~isempty(strfind(mam_list(ii).name, 'R'))
        mam = fliplr(mam);
    end
    
    %segment the breast
    breast_border = segment_breast('image', mam, 'mlo', mlo, 'plot', 0);
    
    %do something with border - (remember to subtract 1 before scaling,
    %then add 1 after)
    
    %set(gcf, 'Name', mam_list(ii).name);
    %saveas(gcf, ['C:\isbe\dev\segmentation\figures\', mam_list(ii).name(1:end-3), 'fig']);
    %close(gcf);
end