function temp_script1(idx)

display(['temp_script1 running: ' datestr(now)]);
%--------------------------------------------------------------------------

iter_per_job = 20;
start_idx = (idx-1)*iter_per_job + 1;
end_idx = min(200, idx*iter_per_job);

for ii = start_idx:end_idx
    %Generate a new random seed (this will be saved)
%     rand_seed = sum(ii*clock);
    load([asymmetryroot 'data/synthetic_data/real512_dt/' zerostr(ii,3) '/rand_seed.mat']);
    rand('twister', rand_seed);
    
    %Make dir to save sampled to
    save_dir = [asymmetryroot 'data/synthetic_data/real512_linop/' zerostr(ii,3) '/'];
    mkdir(save_dir);
    
    %Set up sampling arguments
    display(['generating dataset ' num2str(ii)]);
    sampling_args.num_samples = 2e5;
    sampling_args.bg_dir = [asymmetryroot 'data/synthetic_backgrounds/real512/train/'];
    sampling_args.decomp_type = 'linop';
    sampling_args.num_bgs = 840;
    sampling_args.bg_stem = 'bg';
    sampling_args.bg_zeros = 5;
    sampling_args.detection_type = 'orientation';
    sampling_args.pts_per_image = 500;
    sampling_args.bg_ratio = 1;
    sampling_args.width_range = [4 16];
    sampling_args.orientation_range = [0 360];
    sampling_args.contrast_range = [4 8];
    sampling_args.decay_rate = 4;
    sampling_args.line_type = 'sin';
    sampling_args.normalise = 0;
    sampling_args.num_levels = 5;
    sampling_args.num_angles = 6;
    sampling_args.do_max =  0;
    sampling_args.rotate = 0;
    sampling_args.win_size = 3;
    sampling_args.pca = [];
    sampling_args.use_nag = 0;
    sampling_args.plot = 0;
    sampling_args.save_path = [];
    
    %Sample data
    [X y parameters] = generate_line_training_data(sampling_args); %#ok
    
    %Reshape training data to num_samples*win_size^2*num_bands*num_levels
    X = reshape(X, 2e5, 9, 6, 5); %#ok
    
    %Save the data
    %save([save_dir 'X.mat'], 'X');
    %save([save_dir 'y.mat'], 'y'); 
    save([save_dir 'parameters.mat'], 'parameters');
    save([save_dir 'rand_seed.mat'], 'rand_seed');

    for jj = 1:20
        keep_rows = (1:1e4) + (jj-1)*1e4;
        X_10k = X(keep_rows,:,:,:); %#ok
        save([save_dir 'X_' zerostr(jj,2) '.mat'], 'X_10k');
        clear X_10k;
        
        y_10k = y(keep_rows); %#ok
        save([save_dir 'y_' zerostr(jj,2) '.mat'], 'y_10k');
        clear y_10k;
    end
    clear X y;

end
%--------------------------------------------------------------------------
display(['temp_script1 completed: ' datestr(now)]);

%--------------------------------------------------------------------------
% Old Stuff
%--------------------------------------------------------------------------
% A_mags = zeros(6,6,360);
% A_phases = zeros(6,6,360);
% 
% lev = 3;
% 
% cls_15 = create_gauss_bar(2, 1, 15, 256, 256, 128.5, 128.5);
% dt_15 = dtwavexfm2b(cls_15, 4);
% dts_15 = dt_to_pixel_subset(dt_15, 128.5, 128.5);
% dts_15(imag(dts_15) < 0) = conj(dts_15(imag(dts_15) < 0));
%     
% mags_15 = squeeze(abs(dts_15(:,:,:,lev)));
% phases_15 = squeeze(angle(dts_15(:,:,:,lev)));
% 
% mags_theta = zeros(6,360);    
% phases_theta = zeros(6,360);
% 
% for jj = 1:360
% 
%     cls_theta = create_gauss_bar(2, 1, jj, 256, 256, 128.5, 128.5);
% 
%     dt_theta = dtwavexfm2b(cls_theta, 4);
%     dts_theta = dt_to_pixel_subset(dt_theta, 128.5, 128.5);
%     dts_theta(imag(dts_theta) < 0) = conj(dts_theta(imag(dts_theta) < 0));
% 
%     mags_theta(:,jj) = squeeze(abs(dts_theta(:,:,:,3)));
%     phases_theta(:,jj) = squeeze(angle(dts_theta(:,:,:,3)));
%     
%     A_mags(:,:,jj) = mags_15 / mags_theta(:,jj);
%     A_phases(:,:,jj) = phases_15 / phases_theta(:,jj);
% 
% end
% %save C:\isbe\asymmetry_project\data\misc\A.mat A* 
% 
% figure('name', 'DT magnitudes at centre of rotation');
% for ii = 1:6
%     subplot(2,3,ii);
%     plot(1:360, mags_theta(ii,:)); axis([1 360 0 max(mags_theta(:))]);
%     xlabel('Angle (degrees)');
%     ylabel('DT-CWT magnitude');
%     title(['Band ' num2str(ii)]);  
% end
% 
% figure('name', 'DT phases at centre of rotation');
% for ii = 1:6
%     subplot(2,3,ii);
%     plot(1:360, phases_theta(ii,:)); axis([1 360 -pi pi]);
%     xlabel('Angle (degrees)');
%     ylabel('DT-CWT phase');
%     title(['Band ' num2str(ii)]);
% end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% con_list = dir([asymmetryroot 'data/mammograms/2004_screening/abnormals/mat/*.mat']);
% 
% for ii = 1:length(con_list);
%     m_name = con_list(ii).name(5:10);
%     for g_width = [8 16 32 64]
%         f_name = [asymmetryroot 'data/radial_maps/2004_screening/abnormals/' m_name '_rad_map_' zerostr(g_width,3) '.mat'];
%         
%         if ~exist(f_name, 'file')
%             line_prob = u_load([asymmetryroot 'data/line_maps/2004_screening/abnormals/' m_name '_data.mat']);
%             line_ori = u_load([asymmetryroot 'data/orientation_maps/2004_screening/abnormals/' m_name '_data.mat']);
%     
%             line_ori = imresize(line_ori, 0.5, 'bilinear');
%             line_prob = imresize(line_prob, 0.5, 'bilinear');
%             [angle_bands dist_sum] = radial_line_projection(line_prob, line_ori, [36 1], fspecial('gaussian', [1 5*g_width], g_width)); %#ok
%             save(f_name, 'dist_sum');
%             
%             clear line* angle_bands dist_sum
%             %display(f_name);
%         end
%     end
% end