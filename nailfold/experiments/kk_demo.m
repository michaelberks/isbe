%Script for extracting the raw frames captured during KK Tech's demo of
%their handheld capillaroscopy device

%% 1) Use KK's software to export the movie from it's proprietary format to
%an uncompressed AVI
%-------------------------------------------------------------------------

%% 2) Load in the avi and extract the raw frames
nailfold_dir = 'C:\isbe\nailfold\images\kk_demo\'; %'M:\nailfold\images\kk_demo\';
nailfold_kk = aviread([nailfold_dir 'nailfolds\n3.avi']);
num_frames = length(nailfold_kk);

%Create a holder for the frames - the original frames are 480*752, but
%we're going to discard the first 20 rows because these contain a KK demo
%logo
nailfold_frames = zeros(460, 752, num_frames, 'uint8');

%Make a directory to store the raw frames
mkdir([nailfold_dir 'nailfolds\frames\']);

for ii = 1:num_frames
    %Get frame from the movie structure
    frame = nailfold_kk(ii).cdata(21:end, :);
    
    %Add the frame to our holder
    nailfold_frames(:,:,ii) = frame;
    
    %Save the frame as its own image
    imwrite(frame, [nailfold_dir 'nailfolds\frames\n3_f' zerostr(ii,3) '.bmp']);
end

%save frames in single file for ease of use later
save([nailfold_dir 'nailfolds\n3.mat'], 'nailfold_frames');
%--------------------------------------------------------------------------

%% 3) Try and register the frame into a mosaic using our own algorithm
nailfold_frames = u_load([nailfold_dir 'nailfolds\n3.mat']);

%Test and debug
register_tiles(...
    nailfold_frames(:,:,1:10),...
    'offset_lim', [-100 100; -100 100],...
    'theta_lim', 0,...
    'weights', 'uniform',...
    'debug', 1);
%Run for real
[nailfold_mosaic mosaic_weights offsets thetas] = register_tiles(...
    nailfold_frames,...
    'offset_lim', [-200 200; -100 100],...
    'theta_lim', 0,...
    'weights', 'uniform',...
    'debug', 0);

%Save the mean composition of the registered frames
save([nailfold_dir 'nailfolds\n3_mb.mat'], 'nailfold_mosaic', 'offsets', 'thetas');
imwrite(uint8(nailfold_mosaic), [nailfold_dir 'nailfolds\n3_mb.bmp']);
%--------------------------------------------------------------------------
%% 4) Now try computing the SD of frames across the mosaic rather than
% simply taking the mean
[nailfold_ssq ssq_count] = compute_mosaic_ssq(nailfold_frames, offsets, thetas);

nailfold_sd = sqrt(nailfold_ssq ./ ssq_count - nailfold_mosaic.^2);
nailfold_sd_cap = nailfold_sd;
nailfold_sd_cap(nailfold_sd > 20) = 20;

nailfold_mosaic_sd = nailfold_mosaic - nailfold_sd_cap;
figure; imagesc(nailfold_mosaic_sd); axis image; colormap(gray(256)); caxis([100 240]);

save([nailfold_dir 'nailfolds\n3_std.mat'], 'nailfold_sd', 'nailfold_ssq', 'ssq_count');
imwrite(uint8(nailfold_mosaic_sd), [nailfold_dir 'nailfolds\n3_sd_mb.bmp']);
%--------------------------------------------------------------------------
%% 5) Try and match up the final mosiac to existing 1 of AM
n_am1 = imread('C:\isbe\nailfold\images\ncm\andrea Murray.bmp');
n_am1 = imresize(n_am1(:,:,1), 0.5, 'bilinear');

load([nailfold_dir 'nailfolds\n3_mb.mat'], 'nailfold_mosaic');
n_am2 = uint8(imresize(nailfold_mosaic, 0.5*0.75, 'bilinear')); %New system is x400 magnification, old is x300
clear nailfold_mosaic

%Put them in a common container
sz = max(size(n_am1), size(n_am2));
n_am = zeros(sz(1), sz(2), 2, 'uint8');
n_am(1:size(n_am1,1), 1:size(n_am1,2), 2) = n_am1;
n_am(1:size(n_am2,1), 1:size(n_am2,2), 1) = n_am2;
clear n_am1 n_am2; pack;

masks = n_am > 10 & n_am < 240;
masks(:,:,1) = imerode(masks(:,:,1), strel('disk', 10));
masks(:,:,2) = imerode(masks(:,:,2), strel('disk', 10));

figure; imagesc(masks(:,:,1)); axis image;
figure; imagesc(masks(:,:,2)); axis image;

[nailfold_mosaic mosaic_weights offset theta] = register_tiles(...
    n_am,...
    'tile_masks', masks,...
    'offset_lim', [700 800; -200 200],...
    'theta_range', 15,...
    'max_pts', inf,...
    'weights', 'uniform',...
    'debug', 1);
    

