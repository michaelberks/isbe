clc;
clear;
% close all;
timebar('closeall');

do_profile = false;

if do_profile
    profile clear;
    profile on;
end

% Read input
imgroot = 'U:\projects\nailfold\capture\';
% imgpath = fullfile(imgroot,'2012_10_18\Left\Digit4\x300\seq3');
% imgpath = fullfile(imgroot,'2012_10_22\Left\Digit4\x300\seq2');
% imgpath = fullfile(imgroot,'2013_02_21\Left\Digit4\x300\10_36_55\');
% imgpath = fullfile(imgroot,'2013_06_18\Left.Digit4.x300\10_28_35\');
imgpath = fullfile(imgroot,'2013_07_31\Left.Digit4.x300\12_01_23 (flow 90fps)');

% imgroot = 'N:\Wellcome\videos';
% imgpath = fullfile(imgroot,'short\f72\Left\Digit4\x300');
% imgpath = fullfile(imgroot,'full\f102_v2\Left\Digit4\x300');

use_corrected = true;
if use_corrected
    imgpath = fullfile(imgpath, 'corrected');
end

d = dir(fullfile(imgpath,'frame_*.png'));
% d = d(1:10:end);

n_frames = length(d);
if (n_frames == 0), return; end

frame1 = imread(fullfile(imgpath,d(1).name));
frame1 = mean(frame1,3);

frames = zeros([size(frame1), n_frames], 'uint8');
frames(:,:,1) = uint8(frame1);

for i = 2:n_frames
    frame = imread(fullfile(imgpath,d(i).name));
    frames(:,:,i) = uint8(mean(frame,3));
end

% Register the tiles together (i.e. compute the transformation from one to
% the next over the sequence).
if (1)
    [compound_transforms] = ...
        register_tiles_features(frames, ...
                                'theta_range', [0], ...
                                'offset_lim', 40, ...
                                'debug', false);

    [nailfold_mosaic, mosaic_weights, compound_transforms] = ...
        create_mosaic(frames, compound_transforms);

    regpath = fullfile(imgpath,'registered_g1d');
else
    [nailfold_mosaic mosaic_weights compound_transforms] = ...
        register_tiles_flow(frames);

    regpath = fullfile(imgpath,'registered_flow');
end

if ~exist(regpath,'dir')
    mkdir(regpath);
end
imwrite(uint8(nailfold_mosaic), fullfile(regpath,'mosaic.png'));

% Create batch file that masks the images and creates a movie
maskedpath = fullfile(regpath,'masked');
if ~exist(maskedpath,'dir')
    mkdir(maskedpath);
end

fid = fopen(fullfile(regpath,'make_movie.bat'),'w');
    fprintf(fid,'@echo off\r\n');
    fprintf(fid,'mkdir masked\r\n');
    fprintf(fid,'del masked\\*.png\r\n');
    fprintf(fid,'copy frame*.png masked\r\n');
    fprintf(fid,'mogrify -clip-mask seq_mask.png -threshold 101%%%% masked/frame*.png\r\n');
    fprintf(fid,'ffmpeg -b 1200k -r 30 -i masked/frame_%%%%04d.png -an masked/_movie.mpg\r\n');
fclose(fid);

% Save the images and store the average difference image between every
% frame and the mosaic.
diff_img = write_trans_tiles(frames, compound_transforms, ...
                             regpath, 'frame_', [0 255], nailfold_mosaic);

%% Remove any camera-based artefacts (e.g. dirt on the lens, vignetting).

if ~use_corrected
    figure(2); colormap(gray(256));
        imagesc(diff_img); axis('image');

    imwrite(uint8(255 * normim(diff_img)), fullfile(regpath,'artefacts.png'));

    regpath = fullfile(imgpath,'corrected');
    if ~exist(regpath,'dir')
        mkdir(regpath);
    end
    for f = 1:n_frames
        frame = double(frames(:,:,f)) - diff_img;
        filename = sprintf('frame_%04d.png', f);
        imwrite(uint8(frame), fullfile(regpath, filename));
    end
end

if do_profile
    profile off;
    profile report;
end
