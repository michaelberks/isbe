%----------------------------------------------------------------------
% Efficiency testing:
%----------------------------------------------------------------------
warning('off', 'ASYM:unexpectedArgument');
d_args{1}.decomp_type = {'dt'};
d_args{1}.levels = 1:5;
d_args{1}.feature_shape = 'rect';
d_args{1}.feature_type = 'complex';
d_args{1}.do_max = 0;
d_args{1}.rotate = 0;
d_args{1}.use_nag = 0;
d_args{1}.unwrap_phase = 0;
d_args{1}.interp_mag_phase = 0;     

d_args{2}.decomp_type = {'h2da'};
d_args{2}.sigma_range = [1 2 4 8 16];
d_args{2}.num_angles = 6;
d_args{2}.do_max = 0;
d_args{2}.rotate = 0;

d_args{3}.decomp_type = {'gabor'};
d_args{3}.num_angles = 1;
d_args{3}.sigma_range = [1 2 4 8 16];	
d_args{3}.do_max = 0;
d_args{3}.rotate = 0;
d_args{3}.feature_type = 'complex';

d_args{4}.decomp_type = {'mono'};
d_args{4}.num_levels = 5;
d_args{4}.min_wavelength = 4;
d_args{4}.onf = 0.65;
  
d_args{5}.decomp_type = {'gabori'};
d_args{5}.sigma_range = [1 6];
d_args{5}.num_angles = 5;
d_args{5}.do_max = 0;
d_args{5}.feature_type = 'complex';

d_args{6}.decomp_type = {'g2d'};
d_args{6}.sigma_range = [1 2 4 8 16];
            
d_args{7}.decomp_type = {'g2di'};
d_args{7}.sigma_range = [1 5];
            
d_args{8}.decomp_type = {'h2d'};
d_args{8}.sigma_range = [1 2 4 8 16];

d_args{9}.decomp_type = {'g2da'};
d_args{9}.sigma_range = [1 2 4 8 16];
d_args{9}.num_angles = 6;
d_args{9}.do_max = 0;
d_args{9}.rotate = 0;

num_decomps = 9;

%set common parameters for all decomp types and then compute vector sizes
for i_d = 1:num_decomps
    d_args{i_d}.win_size = 3;
    d_args{i_d}.normalise = 0;
    d_args{i_d}.pca = [];
end
%%

pts_per_im = 1000;
num_ims = 20;
times = zeros(3,num_decomps,2,num_ims);
%%
for i_sz = 1%:3
    sz = 32*(2^i_sz);
        
    for i_im = 1:num_ims    
        im = rand(sz);
        
        %rows = ceil((sz-2)*rand(pts_per_im,1))+1;
        %cols = ceil((sz-2)*rand(pts_per_im,1))+1;
        
        [cols rows] = meshgrid(2:sz-1, 2:sz-1);
        
        for i_decomp = 1%1:num_decomps
            
            tic;
            im_responses = compute_filter_responses(im, d_args{i_decomp}); 
            times(i_sz, i_decomp, 1, i_im) = toc;
            tic;
            sample_image_features(im_responses, rows(:), cols(:), d_args{i_decomp});
            times(i_sz, i_decomp, 2, i_im) = toc;
        end
    end
end
%%
save C:\isbe\asymmetry_project\experiments\times_3.mat times;
%%
figure;
plot([1 2 4 8], mean(times(:,:,1,:), 4));
figure;
plot([1 2 4 8], mean(times(:,:,2,:), 4));

%%
times_dt_f = zeros(100,9);
for i_sz = 0:100
    sz = 60 + 4*i_sz;
    im = rand(sz);
    for i_decomp = 1:num_decomps
        t = zeros(5,1);
        for i_rpt = 1:5
            tic;
            im_responses = compute_filter_responses(im, d_args{i_decomp}); 
            t(i_rpt) = toc;
        end
        if i_sz; times_dt_f(i_sz, i_decomp) = median(t); end
    end
end
%%
m_decomps = zeros(num_decomps,1);

for i_decomp = 1:num_decomps
    x = (60 + 4*(1:100)').^2;
    y = times_dt_f(:,i_decomp);
    p = polyfit(x(26:end), y(26:end), 1);
    m_decomps(i_decomp) = p(1);
    
    ys = polyval(p, x(1));
    ye = polyval(p, x(end));
    
    figure; plot(x, y, 'rx'); hold on;
    plot([x(1) x(end)], [ys ye]);
    
    title(d_args{i_decomp}.decomp_type{1});
end
%%
warning('off', 'ASYM:unexpectedArgument');
times_dt_s = zeros(100,9);
sz = 128;
im = rand(sz);
%%
for i_decomp = [1 5]%1:num_decomps
    im_responses = compute_filter_responses(im, d_args{i_decomp});
    for i_pts = 1:100
        num_pts = 100*i_pts;
        rows = ceil((sz-2)*rand(num_pts,1))+1;
        cols = ceil((sz-2)*rand(num_pts,1))+1;   
    
        t = zeros(5,1);
        for i_rpt = 1:5
            tic;
            sample_image_features(im_responses, rows(:), cols(:), d_args{i_decomp}); 
            t(i_rpt) = toc;
        end
        times_dt_s(i_pts, i_decomp) = median(t);
    end
end
%%
m_decomps_s = zeros(num_decomps,1);
%
for i_decomp = [1 5]%1:num_decomps
    x = (60 + 4*(1:100)').^2;
    y = times_dt_s(:,i_decomp);
    p = polyfit(x(72:end), y(72:end), 1);
    m_decomps_s(i_decomp) = p(1);
    
    ys = polyval(p, x(1));
    ye = polyval(p, x(end));
    
    figure; plot(x, y, 'rx'); hold on;
    plot([x(1) x(end)], [ys ye]);
    
    title(d_args{i_decomp}.decomp_type{1});
end