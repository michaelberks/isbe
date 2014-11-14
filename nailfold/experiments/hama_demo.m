%Script for extracting the raw frames captured during Hamamatsu's deom of
%the Orca R2 and Flash 2.8 cameras

cam_name{1} = 'ORCA Flash 2.8 22ms Exp';
cam_name{2} = 'A Binning 1X1, fast scan';
cam_name{3} = 'A Binning 2X2, fast scan';
cam_name{4} = 'A Binning 1X1, normal scan';
cam_name{5} = 'A Binning 2X2, normal scan';

cam_dir{1} = 'C:\isbe\nailfold\images\camera_demos\oRCA Flash 2.8 Test\';
cam_dir{2} = 'C:\isbe\nailfold\images\camera_demos\ORCA-R2 Test\';
cam_dir{3} = 'C:\isbe\nailfold\images\camera_demos\ORCA-R2 Test\';
cam_dir{4} = 'C:\isbe\nailfold\images\camera_demos\ORCA-R2 Test\';
cam_dir{5} = 'C:\isbe\nailfold\images\camera_demos\ORCA-R2 Test\';

fps = [10 16.2 16.2 9 9];

%% 1) Load in the frames from the multi-tiff files and make
for jj = 1:5
      
    vid_info = imfinfo([cam_dir{jj} cam_name{jj} '.tif']);
    
    counts = zeros(1,4096);
    for ii = 1:length(vid_info);
        flash = imread([cam_dir{jj} cam_name{jj} '.tif'], 'index', ii);
        counts = counts + hist(flash(:), 0:4095);
    end
    figure; bar(0:4095, counts);
    min_g = find(counts, 1, 'first');
    max_g = find(counts, 1, 'last');
    
    if ~isdir([cam_dir{jj} 'frames\']); mkdir([cam_dir{jj} 'frames\']); end
    for ii = 1:length(vid_info);
        vid_frame = double(imread([cam_dir{jj} cam_name{jj} '.tif'], 'index', ii));
        small_frame = uint8(255*(vid_frame - min_g) /(max_g - min_g));
        imwrite(small_frame, [cam_dir{jj} 'frames\' cam_name{jj} sprintf('%03d.png',ii)]);
    end
    
    cmd = ['ffmpeg -f image2 -i "' cam_dir{jj} 'frames\' cam_name{jj} '%03d.png" -r ' sprintf('%3.1f', fps(jj)) ' "' cam_dir{jj} cam_name{jj} '.mp4"'];
    [status result] = system(cmd);
    
end

%% 2) Try and register the frame into a mosaic using our own algorithm
for jj = 2:4
      
    vid_info = imfinfo([cam_dir{jj} cam_name{jj} '.tif']);
    num_frames = length(vid_info);
    nailfold_frames = zeros(vid_info(1).Height, vid_info(1).Width, num_frames, 'uint8');
    
    for ii = 1:num_frames
        nailfold_frames(:,:,ii) = imread([cam_dir{jj} 'frames\' cam_name{jj} sprintf('%03d.png',ii)]);
    end

    %Test and debug
    if 0
        register_tiles(...
            nailfold_frames(:,:,1:10),...
            'offset_lim', [-100 100; -100 100],...
            'theta_lim', 0,...
            'weights', 'uniform',...
            'debug', 1);
    end
    
    %Run for real
    if 1
        [nailfold_mosaic mosaic_weights offsets thetas] = register_tiles(...
            nailfold_frames,...
            'offset_lim', [-200 200; -100 100],...
            'theta_lim', 0,...
            'weights', 'uniform',...
            'debug', 0);
    end

    mkdir([cam_dir{jj} 'mosaics\']);
    %Save the mean composition of the registered frames
    save([cam_dir{jj} 'mosaics\' cam_name{jj} '_mosaic.mat'], 'nailfold_mosaic', 'offsets', 'thetas');
    imwrite(uint8(nailfold_mosaic), [cam_dir{jj} cam_name{jj} '_mosaic.bmp']);
    
    %%Now try computing the SD of frames across the mosaic rather than
    % simply taking the mean
    [nailfold_ssq ssq_count] = compute_mosaic_ssq(nailfold_frames, offsets, thetas);

    nailfold_sd = sqrt(nailfold_ssq ./ ssq_count - nailfold_mosaic.^2);
    nailfold_sd_cap = nailfold_sd;
    nailfold_sd_cap(nailfold_sd > 20) = 20;

    nailfold_mosaic_sd = nailfold_mosaic - nailfold_sd_cap;
    figure; imagesc(nailfold_mosaic_sd); axis image; colormap(gray(256)); caxis([100 240]);

    save([cam_dir{jj} 'mosaics\' cam_name{jj} '_mosaic_sd.mat'], 'nailfold_sd', 'nailfold_ssq', 'ssq_count');
    imwrite(uint8(nailfold_mosaic_sd), [cam_dir{jj} 'mosaics\' cam_name{jj} '_mosaic_sd.bmp']);
end
%%
for jj = 2:5
      
    vid_info = imfinfo([cam_dir{jj} cam_name{jj} '.tif']);
    num_frames = length(vid_info);
    nailfold_frames = zeros(vid_info(1).Height, vid_info(1).Width, num_frames, 'uint8');
    
    for ii = 1:num_frames
        nailfold_frames(:,:,ii) = imread([cam_dir{jj} 'frames\' cam_name{jj} sprintf('%03d.png',ii)]);
    end

    %Save the mean composition of the registered frames
    load([cam_dir{jj} 'mosaics\' cam_name{jj} '_mosaic.mat'], 'offsets', 'thetas');
    
    g_lims = double([min(nailfold_frames(:)) max(nailfold_frames(:))]);
    write_trans_tiles(nailfold_frames, offsets, thetas, [cam_dir{jj} 'reg_frames\'], cam_name{jj}, g_lims);
    
    cmd = ['ffmpeg -f image2 -i "' cam_dir{jj} 'reg_frames\' cam_name{jj} '%03d.png" -r ' sprintf('%3.1f', fps(jj)) ' "' cam_dir{jj} cam_name{jj} '_reg.mp4"'];
    %[status result] = system(cmd);
    display(cmd);
end
%%

%%
jj = 5;
    
load([cam_dir{jj} 'mosaics\' cam_name{jj} '_mosaic_sd.mat'], 'nailfold_sd', 'nailfold_ssq', 'ssq_count');
load([cam_dir{jj} 'mosaics\' cam_name{jj} '_mosaic.mat'], 'nailfold_mosaic', 'offsets', 'thetas');

nailfold_sd_cap = nailfold_sd;
nailfold_sd_cap(nailfold_sd > 20) = 20;
figure; imagesc(nailfold_sd_cap); axis image;
figure; imagesc(nailfold_mosaic); axis image; axis image; colormap(gray(256)); caxis([50 240]);
figure; imagesc(nailfold_mosaic - nailfold_sd_cap); axis image; colormap(gray(256)); caxis([50 240]);

%%
jj = 3;
vid_info = imfinfo([cam_dir{jj} cam_name{jj} '.tif']);
num_frames = 16;%length(vid_info);
rows = vid_info(1).Height / 2;
cols = vid_info(1).Width;
nailfold_frames = zeros(rows, cols, num_frames, 'uint8');

for ii = 1:num_frames
     full_frame = imread([cam_dir{jj} 'frames\' cam_name{jj} sprintf('%03d.png',ii)]);
     nailfold_frames(:,:,ii) = full_frame(rows+1:end,:);
end

dx(:,:,1) = [1 -1; 1 -1] / 4;
dx(:,:,2) = [1 -1; 1 -1] / 4;
dy(:,:,1) = [-1 -1; 1 1] / 4;
dy(:,:,2) = [-1 -1; 1 1] / 4;
dt(:,:,1) = [1 1; 1 1] / 4;
dt(:,:,2) = [-1 -1; -1 -1] / 4;
%%
dx = zeros(41,41,41);
[g,dg] = gaussian_filters_1d(4,20);
gc = g'*dg;
for ii = 1:41
    dx(:,:,ii) = gc * g(ii);
end
dy = permute(dx, [2 1 3]);
dt = permute(dx, [3 2 1]);


laplace = [1 2 1; 2 12 2; 1 2 1]/12;

Ex = convn(nailfold_frames, dx, 'same');
Ey = convn(nailfold_frames, dy, 'same');
Et = convn(nailfold_frames, dt, 'same');
alpha2 = 2^2;

%%
u_flow = zeros(rows, cols);
v_flow = zeros(rows, cols);

for tt = 1:16
    Ex2 = Ex(:,:,tt).^2;
    Ey2 = Ey(:,:,tt).^2;
    Et2 = Et(:,:,tt).^2;

    for itr = 1:32

        u_bar = conv2(u_flow, laplace, 'same');
        v_bar = conv2(v_flow, laplace, 'same');

        aux = (Ex(:,:,tt).*u_bar + Ey(:,:,tt).*v_bar + Et(:,:,tt)) ./ ...
            (alpha2 + Ex2 + Ey2);

        u_flow = u_bar - Ex(:,:,tt) .* aux;
        v_flow = v_bar - Ey(:,:,tt) .* aux;
    end
    
    figure; 
    subplot(2,1,1); quiver(u_flow, v_flow); axis image;
    subplot(2,1,2); imagesc(u_flow.^2 + v_flow.^2); axis image; colorbar;
end
%%
for jj = 1
      
    vid_info = imfinfo([cam_dir{jj} cam_name{jj} '.tif']);
    num_frames = length(vid_info);
    nailfold_frames = zeros(vid_info(1).Height/2, vid_info(1).Width/2, num_frames, 'uint8');
    
    for ii = 1:num_frames
        nailfold_frames(:,:,ii) = imresize(imread([cam_dir{jj} 'frames\' cam_name{jj} sprintf('%03d.png',ii)]), 0.5, 'bilinear');
    end

    %Test and debug
    if 0
        register_tiles(...
            nailfold_frames(:,:,1:10),...
            'offset_lim', [-100 100; -100 100],...
            'theta_lim', 0,...
            'weights', 'uniform',...
            'debug', 1);
        continue;
    end
    
    %Run for real
    [nailfold_mosaic mosaic_weights offsets thetas] = register_tiles(...
        nailfold_frames,...
        'offset_lim', [-200 200; -200 200],...
        'theta_lim', 0,...
        'weights', 'uniform',...
        'debug', 0);

    mkdir([cam_dir{jj} 'mosaics\']);
    %Save the mean composition of the registered frames
    save([cam_dir{jj} 'mosaics\' cam_name{jj} '_mosaic.mat'], 'nailfold_mosaic', 'offsets', 'thetas');
    imwrite(uint8(nailfold_mosaic), [cam_dir{jj} cam_name{jj} '_mosaic.bmp']);
    
    %%Now try computing the SD of frames across the mosaic rather than
    % simply taking the mean
    [nailfold_ssq ssq_count] = compute_mosaic_ssq(nailfold_frames, offsets, thetas);

    nailfold_sd = sqrt(nailfold_ssq ./ ssq_count - nailfold_mosaic.^2);
    nailfold_sd_cap = nailfold_sd;
    nailfold_sd_cap(nailfold_sd > 20) = 20;

    nailfold_mosaic_sd = nailfold_mosaic - nailfold_sd_cap;
    figure; imagesc(nailfold_mosaic_sd); axis image; colormap(gray(256)); caxis([100 240]);

    save([cam_dir{jj} 'mosaics\' cam_name{jj} '_mosaic_sd.mat'], 'nailfold_sd', 'nailfold_ssq', 'ssq_count');
    imwrite(uint8(nailfold_mosaic_sd), [cam_dir{jj} 'mosaics\' cam_name{jj} '_mosaic_sd.bmp']);
    
    g_lims = double([min(nailfold_frames(:)) max(nailfold_frames(:))]);
    write_trans_tiles(nailfold_frames, offsets, thetas, [cam_dir{jj} 'reg_frames\'], cam_name{jj}, g_lims);
    
    cmd = ['ffmpeg -f image2 -i "' cam_dir{jj} 'reg_frames\' cam_name{jj} '%03d.png" -r ' sprintf('%3.1f', fps(jj)) ' "' cam_dir{jj} cam_name{jj} '_reg.mp4"'];
    %[status result] = system(cmd);
    display(cmd);
end
    
    


