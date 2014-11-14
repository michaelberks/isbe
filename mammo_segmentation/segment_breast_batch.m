function [failed_list] = segment_breast_batch(mam_dir, segmentation_dir, image_fmt, mam_list, if_plot, manual)
%SEGMENT_BREAST_BATCH Applies breast segementation to algorithm to all
%mammograms in a directory
%   [failed_list] = segment_breast_batch(mammo_dir, segmentation_dir)
%
% Inputs:
%      mam_dir- Directory containing mammograms to segment
%
%      segmentation_dir- Directory in which to store segmented outputs
%
% Optional inputs (default):
%
%       image_fmt ("tif")- str containing extension type of mammogram images
%
%       mam_list ([])- struct containing specific list of mammograms to
%       process, otherwise all in floder processed
%
%       if_plot (0)- if 1, display debugging visual output from segmentation
%
%       manual (0)- if 1, override the automatic orientation and manually
%       orient the mammogram
%
%
% Outputs:
%      failed_list- List of mammograms for which an error code was returned
%      during segementation
%
%
% Example:
%
% Notes: Assumes mammograms are named containing either their viewcode
% either: LML, LCC, RML, RCC. Will attempt to correctly orient the
% mammogram - right breasts will then be flipped in SEGMENT_BREAST. Images
% are resized to have 1024 rows and an 8-bit gray scale
%
% See also: SEGMENT_BREAST
%
% Created: 10-Jun-2010
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 5125 
% Copyright: (C) University of Manchester

%Check directories are either empty (i.e. use current directory) or filesep
%terminated
if ~isempty(mam_dir) && ~strcmp(mam_dir(end), filesep)
    mam_dir = [mam_dir filesep];
end
if ~isempty(segmentation_dir) && ~strcmp(segmentation_dir(end), filesep)
    segmentation_dir = [segmentation_dir filesep];
end

%Check whether we've been given a specific list of mammograms to process
if ~exist('mam_list', 'var') || isempty(mam_list)
    
    %Check whether a particular image format is specified - default 'tif'
    if ~exist('image_fmt', 'var') || isempty(image_fmt)
        image_fmt = 'tif';
    end
    
    mam_list_s = dir([mam_dir '*' image_fmt]);
    num_mammos = length(mam_list_s);
    mam_list = cell(num_mammos,1);
    for i_m = 1:num_mammos
        mam_list{i_m} = mam_list_s(i_m).name;
    end
else
    num_mammos = length(mam_list);
end

%Check whether we've been asked to give visual output
if nargin < 5
    if_plot = 0;
end

%Check whether we've been asked to give visual output
if nargin < 6
    manual = 0;
end

%Check whether segmentation directory exists
if ~isdir(segmentation_dir)
    mkdir(segmentation_dir);
end    

%pre-allocate failures list
failed_list = [];
for ii = 1:num_mammos
    try
        %load in the mammogram and resize to have 1024 rows
        mam = imresize(imread([mam_dir mam_list{ii}]), [1024 NaN], 'bilinear');
        
        %Check mammogram is portrait orientated
        [r c] = size(mam);
        if c > r
            mam = rot90(mam);
        end
        
        % resize to have 1024 rows
        mam = imresize(mam, [1024 NaN], 'bilinear');
        
        %Check mammo format - if uint16 need to divide by 256
        if isa(mam, 'uint16')
            mam = mam / 256;
        end
        
        %Check whether mammogram is MLO or CC and right or left
        mlo = ~isempty(strfind(mam_list{ii}, 'LML')) || ~isempty(strfind(mam_list{ii}, 'RML'));
        right = ~isempty(strfind(mam_list{ii}, 'RML')) || ~isempty(strfind(mam_list{ii}, 'RCC'));
        
        %Check mammogram is correctly oriented
        if manual
            f1 = figure; imagesc(mam); axis image; colormap(gray(256));
            title(mam_list{ii});
            button = questdlg('Is the mammogram correctly oriented?','Manual orientation','Yes','No','Yes');
            close(f1);
            if strcmpi(button, 'no')
                mam = rot90(mam, 2);
            end

        else
            mam = orientate_mam(mam, right);
        end
        
        %Complete the segmentation
        [segmentation] = segment_breast('mammo', mam, 'mlo', mlo, 'right', right, 'plot', if_plot); %#ok
        
        %Save to the output directory
        save([segmentation_dir mam_list{ii}(1:end-4) '_segmentation.mat'], 'segmentation');
    catch err
        %add failed mammogram name to failures list
        failed_list(end+1,1).name = mam_list{ii}; %#ok;
        display(['Skipping ', mam_list{ii} '. ' err.message]);
    end
end

%Write out summary of results
num_failed = length(failed_list);
num_passed = length(mam_list) - num_failed;

results_file = [segmentation_dir '\segmentation_results.txt'];

fid1 = fopen(results_file, 'wt');
fprintf(fid1, '%s \n', ['Segmentation performed on ' num2str(num_mammos)]);
fprintf(fid1, '%s \n', [num2str(num_passed) ' processed successfully']);
fprintf(fid1, '%s \n', [num2str(num_failed) ' failed to process']);
fprintf(fid1, '%s \n', 'List of failed mammograms:');

for i_f = 1:num_failed
    fprintf(fid1, '%s \n', failed_list(i_f).name);
end                
fclose(fid1);
                

function mam = orientate_mam(mam, right)

    %Compute the sum of intensities in the 2nd and 3rd quartile of columns
    [r c] = size(mam);
    
    r1 = round(r / 3);
    r2 = 2*r1;
    c_half = round(c / 2);
    
    mam_mean = mean(mam(:));
    centre_mean = mean(mam(r1:r2,:));
    
    sum1 = sum(centre_mean(1:c_half) > mam_mean);
    sum2 = sum(centre_mean(c_half+1:end) > mam_mean);
    
    %if sum1 > sum2 then the breast is on the left of the image
    if (sum1 > sum2) == right
        mam = rot90(mam, 2);
    end
    


