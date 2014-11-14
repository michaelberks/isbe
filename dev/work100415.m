test_dir = 'C:\isbe\dev\classification\data\testimage_contrast0to8_multibars_sin\test_images\';
prob_dir = 'C:\isbe\dev\classification\data\testimage_contrast0to8_multibars_sin\probability_images\';

bins = (0:200)/200;

counts_line = zeros(201,1);
counts_bg = zeros(201,1);

for ii = 1:100
    load([test_dir '\test_image' zerostr(ii,3) '.mat']);
    prob_image = u_load([prob_dir '182270\probability_image' zerostr(ii, 3)]);
    counts_line = counts_line + histc(prob_image(label_centre), bins);
    counts_bg = counts_bg + histc(prob_image(~label_centre), bins);
end

figure; 
subplot(1,2,1); bar(bins, counts_line);
subplot(1,2,2); bar(bins, counts_bg);


%%
figure;
for ii = 1:20
    offset = 80;
    load([test_dir '\test_image' zerostr(ii+offset,3) '.mat']);
    prob_image = u_load([prob_dir '182270\probability_image' zerostr(ii+offset, 3)]);
    im_counts = histc(prob_image(~label_centre), bins);
    subplot(4,5,ii); bar(bins, im_counts); axis([0 1 0 max(im_counts(:))]);
end
%%
%--------------------------------------------------------------------------
%%
prob_dir = 'M:\chen\data\testimage_contrast1to8_exprnd_ellipse\probability_image\';
test_dir = 'M:\chen\data\testimage_contrast1to8_exprnd_ellipse\';

bins = (0:200)/200;

counts_line = zeros(201,1);
counts_bg = zeros(201,1);
for ii = 1:100
    load([test_dir '\label' num2str(ii) '.mat']);
    prob_image = u_load([prob_dir '182393\probability_image' zerostr(ii, 3)]);
    counts_line = counts_line + histc(prob_image(multibar_labelcentre), bins);
    counts_bg = counts_bg + histc(prob_image(~multibar_labelcentre), bins);
end

figure; 
subplot(1,2,1); bar(bins, counts_line);
subplot(1,2,2); bar(bins, counts_bg);
%%

figure;
for ii = 1:20
    offset = 80;
    load([test_dir '\label' num2str(ii+offset) '.mat']);
    prob_image = u_load([prob_dir '182270\probability_image' zerostr(ii+offset, 3)]);
    im_counts = histc(prob_image(multibar_labelcentre), bins);
    subplot(4,5,ii); bar(bins, im_counts); axis([0 1 0 max(im_counts(:))]);
end
%%
bg_dir = 'M:\chen\data\smooth512x512_patches\train\';
bg = u_load([bg_dir 'bg00001.mat']);

for width = 4:16
    [bar, label, label_centre] = create_ellipse_bar(width/2, 4, 24, 512, 512, 256, 256);

    [rr cc] = find(label_centre);
    %[rows cols] = find(label_centre);
    % %Make copies of sample rows and cols at positions of local window patch
    % rr = repmat(rows*ones(1,win_size) + ones(num_samples_image,1)*win_idx, 1, win_size);
    % cc = kron(cols*ones(1,win_size) + ones(num_samples_image,1)*win_idx, ones(1,win_size));

    dt_bg = dtwavexfm2b(bg, 6);

    %Get interpolated dual-tree coefficients
    dt_samples_bg = dt_to_pixel_subset(dt_bg, rr, cc);

    figure; hold all;
    ii = 1;
    for c = [0 0.5 1 2 4 8 12 16]

        [bar, label, label_centre] = create_ellipse_bar(width/2, c, 24, 512, 512, 256, 256);
        dt_bar = dtwavexfm2b(bg+bar, 6);

        %Get interpolated dual-tree coefficients
        dt_samples_bar = dt_to_pixel_subset(dt_bar, rr, cc);

        display(['contrast = ' num2str(c)]);
        %
        mean_diffs = zeros(6,1);
        max_diffs = zeros(6,1);
        for lev = 1:6
            lev_cols = 6*(lev-1)+1:6*lev;
            mean_diffs(lev) = mean(mean(abs(dt_samples_bar(:,lev_cols)-dt_samples_bg(:,lev_cols))));
            max_diffs(lev) = max(max(abs(dt_samples_bar(:,lev_cols)-dt_samples_bg(:,lev_cols))));
            %display(['mean diff for level ' num2str(lev) ' = ' num2str(mean_diffs(lev))]);
            %display(['max diff for level ' num2str(lev) ' = ' num2str(max_diffs(lev))]);
        end
        legend_text{ii} = ['c = ', num2str(c)]; %#ok
        plot(1:6, max_diffs);
        ii = ii+1;
    end
    legend(legend_text);
    title(['Mean change in DT coefficients at centre line - bar width = ', num2str(width)]);
    xlabel('DT level');
    ylabel('Mean change');
end
%%
rand('twister', 5489);
sampling_method_args.bg_dir = 'M:\chen\data\smooth512x512_patches\train\';
sampling_method_args.save_path = [];
sampling_method_args.num_samples = 1e3;
sampling_method_args.width_range = [0 8];
sampling_method_args.contrast_range = [0 8];
sampling_method_args.win_size = 3;
sampling_method_args.num_levels = 5;
sampling_method_args.num_angles = 8;
[X_test3 y_test3] = sample_training_data_linop(sampling_method_args);
%
rand('twister', 5489);
sampling_method_args.bg_dir = '\\isbe-san1\mberks\dev\classification\data\smooth512x512_patches\train\';
sampling_method_args.save_path = [];
sampling_method_args.num_samples = 1e4;
sampling_method_args.width_range = [0 8];
sampling_method_args.contrast_range = [0 8];
sampling_method_args.win_size = 1;
sampling_method_args.num_levels = 5;
sampling_method_args.num_angles = 8;
[X_test1 y_test1] = sample_training_data_linop(sampling_method_args);

%%
X_test1a = X_test3(:,5:9:end);
display(max(abs(X_test1a(:) - X_test1(:))));
display(max(abs(y_test3(:) - y_test1(:))));
%%
mkdir \\isbe-san1\mberks\dev\classification\data\smooth512x512_patches\train\
bg_list = dir('M:\chen\data\smooth512x512_patches\train\*.mat');
r_idx = randsample(1:length(bg_list), 1000);
%%
for ii = 1:1000
    bg = u_load(['M:\chen\data\smooth512x512_patches\train\' bg_list(r_idx(ii)).name]);
    save(['\\isbe-san1\mberks\dev\classification\data\smooth512x512_patches\train\bg' zerostr(ii, 4)], 'bg');
    clear bg;
end
%%
x = 4; mu = 2;
y = expcdf(x, mu)
y2 = 1 - exp(-x / mu)
x2 = -mu*log(1-y);

uni_sample = rand(1e5, 1);
x_sample = -mu*log(1-uni_sample);
x_sample2 = exprnd(mu, 1e5, 1);
x_sample3 = exp_rand_sample(mu, [1e5 1]);
x_sample4 = zeros(1e5,1);
for ii = 1:1e5; x_sample4(ii) = exp_rand_sample(mu); end