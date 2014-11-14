% generating_regions_script

% A script describing how we generated the mammogram regions we've used in the
% project, and also the associated wavelet decompositions

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% NORMAL REGIONS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Select a set of normal regions from mammograms using the function:
if 0;
    mb_sample_valid_region('MammoDir', 'C:\isbe\mammograms\new_CAD\BMP_2004\',...
        'MammoFiles', normal_list, ...
        'SaveDir', 'C:\isbe\dev\background\images\normal1024\',...
        'WindowSize', 1024);
end

%Where normal list a list of file names for the normal disease-free mammograms
%Note this function has a random element, so to obtain the same set of
%regions set the random seed prior to executing

%Having generated 1024x1024 regions at full resolution, we downsample each
%region by 2 to 512x512 and save again

%The code also displays the effect on the local image variance of different
%downsampling interpolation methods - we choose 'bilinear' interpolation

normal_list = dir('C:\isbe\dev\background\images\normal1024\*.bmp');

for ii = 1:length(normal_list);
    
    %load in normal regions
    normal_region_1024 = imread(['C:\isbe\dev\background\images\normal1024\', normal_list(ii).name]);
    
    %Resize by a half using 3 method option
    normal_region_nearest = imresize(normal_region_1024, 0.5, 'nearest');
    normal_region_bilinear = imresize(normal_region_1024, 0.5, 'bilinear');
    normal_region_bicubic = imresize(normal_region_1024, 0.5, 'bicubic');
    
    % Compute and display the local image variance
    [norm_var_1024] = image_stats(normal_region_1024, 5);
    [norm_var_nearest] = image_stats(normal_region_nearest, 5);
    [norm_var_bilinear] = image_stats(normal_region_bilinear, 5);
    [norm_var_bicubic] = image_stats(normal_region_bicubic, 5);
    
    display([norm_var_1024 norm_var_nearest norm_var_bilinear norm_var_bicubic]);
    
    %Save the region for the bilinear downsampling
    imwrite(normal_region_bilinear, ['C:\isbe\dev\background\images\normal512\', normal_list(ii).name]);
    
end
%%
% Now build the dual_tree decompositions for the regions
mb_build_dual_tree(...
    'ImageDir', 'C:\isbe\dev\background\images\normal1024\',...
    'OutputDir', 'C:\isbe\dev\background\dual_tree\normal1024\', 'Levels', 7);

% Now build the dual_tree decompositions for the regions
mb_build_dual_tree(...
    'ImageDir', 'C:\isbe\dev\background\images\normal512\',...
    'OutputDir', 'C:\isbe\dev\background\dual_tree\normal512\', 'Levels', 6);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% MASS REGIONS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The mass regions are created by subtracting the predicted mass values
% from the original mass region. Note so far we have created these for the
% 'good' mass background
% Again, we also downsample each region by a factor two

load('C:\isbe\dev\files\bg_files.mat');
%mkdir C:\isbe\dev\background\images\mass_2
%mkdir C:\isbe\dev\background\images\mass
for ii = 1:length(bg_files);
    
    load(['C:\isbe\dev\masses\', bg_files(ii).name]);
    mass_bg = double(mass.mass_ROI) - mass.subtract_ROI;
    
    %Reduce resolution by a half
    mass_bg2 = imresize(mass_bg, 0.5, 'bilinear');    
    imwrite(uint8(mass_bg), ['C:\isbe\dev\background\images\mass\mass', zerostr(ii,3), '.bmp']);
    imwrite(uint8(mass_bg2), ['C:\isbe\dev\background\images\mass_2\mass', zerostr(ii,3), '.bmp']);
    
end
%
mb_build_dual_tree(...
    'ImageDir', 'C:\isbe\dev\background\images\mass\',...
    'OutputDir', 'C:\isbe\dev\background\dual_tree\mass\', 'Levels', 7);

mb_build_dual_tree(...
    'ImageDir', 'C:\isbe\dev\background\images\mass_2\',...
    'OutputDir', 'C:\isbe\dev\background\dual_tree\mass_2\');
%%