orientation_range = [0 360];
contrast_range = [4 8];
decay_rate = 4;
width_range = [4 16];
%%  
for ii = 1:20
    
    mu = (contrast_range(2) - contrast_range(1)) / (2*log(decay_rate));
    contrast = contrast_range(1) + exp_rand_sample(mu);
    bar = create_sin_curve(3, contrast, 2048, 23, 0, 512, 512, 200, 137);
    figure; imagesc(sampled_window + bar); axis image; colormap(gray(256));
end
%%
for ii = 1:20
    %width = width_range(1) + (width_range(2)-width_range(1))*rand;
    mu = (width_range(2) - width_range(1)) / (2*log(decay_rate));
    width = width_range(1) + exp_rand_sample(mu);
    bar = create_sin_curve(width, 8, 2048, 23, 0, 512, 512, 200, 137);
    figure; imagesc(sampled_window + bar); axis image; colormap(gray(256));
end
%%
for ii = 1:20
    squash = rand;
    bar = create_sin_curve(4, 4, 2048, 23, squash, 512, 512, 200, 137);
    figure; imagesc(sampled_window + bar); axis image; colormap(gray(256));
end
%%
for ii = 1:20
    radius = 2^(8*rand + 8);
    bar = create_sin_curve(4, 4, radius, 23, 0, 512, 512, 200, 137);
    figure; imagesc(sampled_window + bar); axis image; colormap(gray(256));
end
%%
for ii = 1:20
    bg_num = ceil(rand*940);
    bg = u_load(['C:\isbe\asymmetry_project\data\synthetic_backgrounds\real512\bg' zerostr(bg_num, 5) '.mat']);
    
    orientation = orientation_range(1) + (orientation_range(2)-orientation_range(1))*rand;
    
    radius = 2^(8*rand + 8);
    
    mu = (width_range(2) - width_range(1)) / (2*log(decay_rate));
    width = width_range(1) + exp_rand_sample(mu);

    mu = (contrast_range(2) - contrast_range(1)) / (2*log(decay_rate));
    contrast = contrast_range(1) + exp_rand_sample(mu);
    
    squash = rand;
    
    bar = create_sin_curve(width/2, contrast, radius, orientation, squash, 512, 512, 256, 256);
    figure; 
    subplot(1,2,1); imagesc(bg + bar); axis image; colormap(gray(256));
    subplot(1,2,2); imagesc(bar); axis image; colormap(gray(256));
end
%%
mu = (width_range(2) - width_range(1)) / (2*log(decay_rate));
widths = width_range(1) + exp_rand_sample(mu, [1e4, 1]);
figure; hist(widths/2, 100);
%%
mu = (contrast_range(2) - contrast_range(1)) / (2*log(decay_rate));
contrasts = contrast_range(1) + exp_rand_sample(mu, [1e4, 1]);
figure; hist(contrasts, 100);
%%
generate_synthetic_curve_images(...
    'num_images', 100, ...
    'save_dir', 'real512\',... % the mandatory arguments
    'bg_dir', 'C:\isbe\asymmetry_project\data\synthetic_backgrounds\real512\');