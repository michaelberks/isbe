
%mam_dir = 'C:\isbe\density\mammograms\JamieToMike\DigBatch5_anon\';
%mam_list = dir([mam_dir '*.tif']);

mam_dir = 'C:\isbe\density\mammograms\MPhys\';
mam_list = dir([mam_dir 'small\*.mat']);
seg_dir = [mam_dir 'segmentations\'];

%
for ii = 101:114%0%length(mam_list)
    try
        seg = u_load([seg_dir mam_list(ii).name(1:end-4) '_segmentation.mat']);
    catch
        display('segmentation not found');
        continue;
    end
    %mam = imread([mam_dir mam_list(ii).name]);
    %[r c] = size(mam);
    %if c > r
    %    mam = rot90(mam);
    %end
    %mam = imresize(mam, seg.size, 'nearest');
    mam = u_load([mam_dir 'small\' mam_list(ii).name]);
    right = ~isempty(strfind(mam_list(ii).name, 'RML')) || ~isempty(strfind(mam_list(ii).name, 'RCC'));
    mam = orientate_mammogram(mam, right);
    
    figure; imgray(mam);
    
    plot(seg.breast_border(:,1), seg.breast_border(:,2), 'g');
    plot(seg.breast_border(seg.breast_air,1), seg.breast_border(seg.breast_air,2), 'rx');
end
%%
mkdir([mam_dir 'small']);
for ii = 1:length(mam_list)
    
    mam = imread([mam_dir mam_list(ii).name]);
    [r c] = size(mam);
    if c > r
        mam = rot90(mam);
    end
    mam = imresize(mam, [1024 nan], 'nearest');
    save([mam_dir 'small\' mam_list(ii).name(1:end-4) '.mat'], 'mam');
end
%%
mkdir([mam_dir 'jpg']);
for ii = 1:100%0%length(mam_list)
    mam = u_load([mam_dir 'small\' mam_list(ii).name]);
    mam = uint8(double(mam) / 2^8);
    imwrite(mam, [mam_dir 'jpg\' mam_list(ii).name(1:end-4) '.jpg']);
end