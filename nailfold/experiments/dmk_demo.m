%Script for extracting the raw frames captured during Hamamatsu's deom of
%the Orca R2 and Flash 2.8 cameras

cam_dir{1} = 'C:\isbe\nailfold\images\camera_demos\2012_08_06\mike_dmk_camera\enhanced\seq1\';
cam_dir{2} = 'C:\isbe\nailfold\images\camera_demos\2012_08_06\rohma_dmk_camera\enhanced\seq1\';
cam_dir{3} = 'C:\isbe\nailfold\images\camera_demos\2012_08_06\rohma_dmk_camera\enhanced\seq2\';
cam_dir{4} = 'C:\isbe\nailfold\images\camera_demos\2012_08_06\rohma_dmk_camera\enhanced\seq3\';
cam_dir{5} = 'C:\isbe\nailfold\images\camera_demos\2012_08_06\rohma_dmk_camera\enhanced\seq4\';

fps = [1 1 1 1 1] * 60;

%% 1) Try and register the frame into a mosaic using our own algorithm

for jj = 2:5
    
    frame_list = dir([cam_dir{jj} '*.png']);
    
    num_frames = length(frame_list);
    nailfold_frames = zeros(480, 640, num_frames, 'uint8');
    
    for ii = 1:num_frames
        nailfold_frames(:,:,ii) = imread([cam_dir{jj} frame_list(ii).name]);
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
    save([cam_dir{jj} 'mosaics\mosaic.mat'], 'nailfold_mosaic', 'offsets', 'thetas');
    imwrite(uint8(nailfold_mosaic), [cam_dir{jj} 'mosaic.bmp']);
    
    %%Now try computing the SD of frames across the mosaic rather than
    % simply taking the mean
    [nailfold_ssq ssq_count] = compute_mosaic_ssq(nailfold_frames, offsets, thetas);

    nailfold_sd = sqrt(nailfold_ssq ./ ssq_count - nailfold_mosaic.^2);
    nailfold_sd_cap = nailfold_sd;
    nailfold_sd_cap(nailfold_sd > 20) = 20;

    nailfold_mosaic_sd = nailfold_mosaic - nailfold_sd_cap;
    figure; imagesc(nailfold_mosaic_sd); axis image; colormap(gray(256));

    save([cam_dir{jj} 'mosaics\mosaic_sd.mat'], 'nailfold_sd', 'nailfold_ssq', 'ssq_count');
    imwrite(uint8(nailfold_mosaic_sd), [cam_dir{jj} 'mosaics\mosaic_sd.bmp']);
end
%%
for jj = 2:5
      
    frame_list = dir([cam_dir{jj} '*.png']);
    
    num_frames = length(frame_list);
    nailfold_frames = zeros(480, 640, num_frames, 'uint8');
    
    for ii = 1:num_frames
        nailfold_frames(:,:,ii) = imread([cam_dir{jj} frame_list(ii).name]);
    end

    %Save the mean composition of the registered frames
    load([cam_dir{jj} 'mosaics\mosaic.mat'], 'offsets', 'thetas');
    
    g_lims = double([min(nailfold_frames(:)) max(nailfold_frames(:))]);
    mkdir([cam_dir{jj} 'reg_frames\']);
    write_trans_tiles(nailfold_frames, offsets, thetas, [cam_dir{jj} 'reg_frames\'], 'frame_r_', g_lims);
    
    cmd = ['ffmpeg -f image2 -i "' cam_dir{jj} 'reg_frames\frame_r_%03d.png" -r ' sprintf('%3.1f', fps(jj)) ' "' cam_dir{jj} 'frames_reg.mp4"'];
    [status result] = system(cmd);
    %display(cmd);
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
%%
fname{1} = 'C:\isbe\nailfold\images\camera_demos\2012_08_06\mike_mosaic\Visit1\Lhand\Digit4\X300\FullResMosaic\STUDENT2V1LD4X3LrgMosaic';
fname{2} = 'C:\isbe\nailfold\images\camera_demos\2012_08_06\mike_mosaic\Visit1\Rhand\Digit5\X300\FullResMosaic\STUDENT2V1RD5X3LrgMosaic';
fname{3} = 'C:\isbe\nailfold\images\camera_demos\2012_08_06\mike_mosaic\Visit1\Lhand\Digit5\X300\FullResMosaic\STUDENT2V1LD5X3LrgMosaic';
fname{4} = 'C:\isbe\nailfold\images\camera_demos\2012_08_06\rohma_mosaic\Visit1\Lhand\Digit4\X300\FullResMosaic\STUDENTV1LD4X3LrgMosaic';
for ii = 1:length(fname)
    nailfold = imread([fname{ii} '.bmp']);
    nailfold_mask = nailfold < 250;
    clims = prctile(double(nailfold(nailfold_mask)), [1 99]);
    imwrite((double(nailfold)-clims(1)) / diff(clims), [fname{ii} '_e.bmp']);
end


    
    


