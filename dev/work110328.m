
dt_coeffs = zeros(360,24);
sampling_args=[];
sampling_args.feature_shape = 'rect';
sampling_args.feature_type = 'conj';
sampling_args.levels = 4:5;
sampling_args.win_size = 1;

for ori = 1:360
    
    bar = create_gauss_bar(2, 1, ori-0.5, 256, 256, 128, 128);
    dt = compute_dual_tree(bar, 5, 0);
    dt_coeffs(ori,:) = sample_dt_data(dt, 126, 127, sampling_args);
end
%
% x = cos(linspace(0, 2*pi, 360)');
% y = sin(linspace(0, 2*pi, 360)');
x = cosd((1:360)-0.5)';
y = sind((1:360)-0.5)';
figure;
for band = 1:6
    subplot(2,3,band);
    plot(x .* dt_coeffs(:,band+6), y .* dt_coeffs(:,band+6));
    axis([-1 1 -1 1]);
end

figure;
for band = 1:6
    subplot(2,3,band); hold on;
    plot(x .* dt_coeffs(:,band+18), y .* dt_coeffs(:,band+18), 'b.');
    plot(pi*x .* dt_coeffs(:,band+6), pi*y .* dt_coeffs(:,band+6),'r-');
%     plot(y .* dt_coeffs(:,band+18), x .* dt_coeffs(:,band+18),'r-');
	xx = cosd((30*band)-15); yy = sind((30*band)-15);
	plot([xx,-xx],[yy,-yy],'g');
	plot([yy,-yy],[-xx,xx],'g:');
    axis('equal',[-pi pi -pi pi]);
end
%%
[X y] = sample_saved_dt_line_data (...
    'saved_data_dir', 'C:\isbe\asymmetry_project\data\precomputed_data\real512_dt\',...
    'num_samples', 1e5,...
    'task_id', 1,...
	'pts_per_image', [], ...
    'id_offset', 0,...
    'feature_shape', 'rect',...
    'feature_type', 'conj',...
    'num_levels', 5,...
    'do_max', 0,...
    'rotate', 0,...
    'win_size', 1,...
    'pca', []);

%
y_deg = mod(round(90*angle(y)/pi) - 1,180) + 1;

X_means = zeros(180, 60);
X_deg_counts = sparse(y_deg, 1, 1, 180, 1);
for ii = 1:60
    X_deg_sum = sparse(y_deg, 1, X(:,ii), 180, 1);
    X_means(:,ii) = full(X_deg_sum ./ X_deg_counts);
end
clear X;
%%
x = cosd((1:180))';
y = sind((1:180))';
valid_idx = ~isnan(X_means(:,1));
figure;
for band = 1:6
    subplot(2,3,band);
    plot(x(valid_idx) .* X_means(valid_idx,band+18), y(valid_idx) .* X_means(valid_idx,band+18));
    axis([-10 10 -10 10]);
end
%
figure;
for band = 1:6
    subplot(2,3,band);
    plot(x(valid_idx) .* X_means(valid_idx,band+48), y(valid_idx) .* X_means(valid_idx,band+48));
    axis([-pi pi -pi pi]);
end
%
figure;
for band = 1:6
    subplot(2,3,band); hold on;
    plot(x(valid_idx) .* X_means(valid_idx,band+18), y(valid_idx) .* X_means(valid_idx,band+18));
    plot(x(valid_idx) .* X_means(valid_idx,band+48), y(valid_idx) .* X_means(valid_idx,band+48), 'r');
    axis([-10 10 -10 10]);
end

figure;
for band = 1:6
    subplot(2,3,band); hold on;
    plot(x(valid_idx) .* (X_means(valid_idx,band+18) - X_means(valid_idx,band+48)),...
        y(valid_idx) .* (X_means(valid_idx,band+18) - X_means(valid_idx,band+48)));
    
    axis([-10 10 -10 10]);
end
%
figure; hold all;
for band = 1:6
    plot(x(valid_idx) .* X_means(valid_idx,band+48), y(valid_idx) .* X_means(valid_idx,band+48));
end
axis([-pi pi 0 pi]); axis equal;
legend({'band 1', 'band 2', 'band 3', 'band 4', 'band 5', 'band 6'});
%%
ori_list = dir('C:\isbe\asymmetry_project\data\orientation_maps\g2d\2004_screening_processed\normals\*_ori.mat');
for ii = 1:length(ori_list)
    display(['processing map ' num2str(ii)]);
    ori_map = load_uint8(['C:\isbe\asymmetry_project\data\orientation_maps\g2d\2004_screening_processed\normals\' ori_list(ii).name]);
    ori_map = mod(pi*ori_map/180, pi);
    save_uint8(['C:\isbe\asymmetry_project\data\orientation_maps\g2d\2004_screening_processed\normals\' ori_list(ii).name], ori_map);
end



    
    