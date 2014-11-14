%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% This script provides code to segment the mammograms for the Bolton
% density project.
%
% See the notes in segment_breast.m for more details on resizing and
% correctly aligning the mammograms prior to segmentation

%%
%--------------------------------------------------------------------------
% This section performs the segmentation
%--------------------------------------------------------------------------

%Set directory in which mammograms are stored and get list of all images
mam_dir = 'F:\';
mam_list = dir('F:\*.tif');
%Set directory in which to store segmentation results
segment_dir = 'C:\isbe\dev\segmentation\';

%Go through each image
for ii = 1:length(mam_list);
    try
        %Perform segmentation in catch loop so we can batch run without
        %breaking down for errors
        
        %Work out if mammogram is MLO or CC (note we do this from the file
        %name so mammograms *must* be correctly named)
        mlo = ~isempty(strfind(mam_list(ii).name, 'ML'));
        
        %Load the mammogram and resize to that required by segmentation
        %algorithm
        mam = double(imresize(imread([mam_dir, mam_list(ii).name]) / 256, [1024 NaN], 'bilinear'));
        
        %Rotate and flip the mammogram (if necessary) depending on whether
        %the mammogram is a) Left or Right b) 18x24 or 24x30
        % Again, we do this from the name, so mammograms *must* be
        % correctly named
        if isempty(strfind(mam_list(ii).name, '2430'))
            if ~isempty(strfind(mam_list(ii).name, 'R'))
                mam = rot90(mam, 1);
                mam = fliplr(mam);
            else
                mam = rot90(mam, -1);
            end
        else
            if ~isempty(strfind(mam_list(ii).name, 'R'))
                mam = fliplr(mam);
            else
                mam = rot90(mam, 2);
            end
        end
        
        %Peform the segmentation
        segmentation.breast_border = segment_breast('image', mam, 'mlo', mlo, 'plot', 0);
        segmentation.size = size(mam);
        
        %Save the results
        save([segment_dir, mam_list(ii).name(1:end-4), '_segmentation.mat'], 'segmentation');
    catch
        
        %Display any mamograms not succesfully segmented
        display(['Skipping ', mam_list(ii).name]);
    end
end
%%
%--------------------------------------------------------------------------
% This section loads in images and displays the segmented breast border for 
% visual evaluation
%--------------------------------------------------------------------------
%Set directory in which mammograms are stored and get list of all images
mam_dir = 'F:\';
mam_list = dir('F:\*.tif');
%Set directory in which to store segmentation results
segment_dir = 'C:\isbe\dev\segmentation\';

%Go through each image
for ii = 1:length(mam_list);
    try
        %Perform segmentation in catch loop so we can batch run without
        %breaking down for errors
        
        %Work out if mammogram is MLO or CC (note we do this from the file
        %name so mammograms *must* be correctly named)
        mlo = ~isempty(strfind(mam_list(ii).name, 'ML'));
        
        %Load the mammogram and resize to that required by segmentation
        %algorithm
        mam = double(imresize(imread([mam_dir, mam_list(ii).name]) / 256, [1024 NaN], 'bilinear'));
        
        %Rotate and flip the mammogram (if necessary) depending on whether
        %the mammogram is a) Left or Right b) 18x24 or 24x30
        % Again, we do this from the name, so mammograms *must* be
        % correctly named
        if isempty(strfind(mam_list(ii).name, '2430'))
            if ~isempty(strfind(mam_list(ii).name, 'R'))
                mam = rot90(mam, 1);
                mam = fliplr(mam);
            else
                mam = rot90(mam, -1);
            end
        else
            if ~isempty(strfind(mam_list(ii).name, 'R'))
                mam = fliplr(mam);
            else
                mam = rot90(mam, 2);
            end
        end
        
        %Load in the segmentation
        load(['C:\isbe\dev\segmentation\', mam_list(ii).name(1:end-4), '_segmentation.mat'], 'segmentation');
        
        %Display the mammogram and the segmented breast border
        figure; imagesc(mam); axis image; hold on;
        plot(segmentation.breast_border(:,1), segmentation.breast_border(:,2), 'r-');
        plot(segmentation.breast_border(:,1), segmentation.breast_border(:,2), 'rx');
    catch
        %Display any mammograms we skip (probably because there was an
        %error in the segmentation so no file was saved)
        display(['Skipping ', mam_list(ii).name]);
    end
end