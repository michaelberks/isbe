%
% Weights script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Runs code needed to compute the weights used in the combined appearance
% models and the errors associated with each set of weights
% methods   1: mean sd per mode
%           2: mean variance pre mode
%           3: Finite difference method
%           4: Optimisation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% How to compute weights:

%
% Load model
load C:\isbe\dev\mass_model\models\model_i500_50K.mat
%

% 1.
% existing model weights are mean sd - normalise by scale weight
%cd C:\isbe\dev\weights
weights_sd = [mass_model.W_shape mass_model.W_tex] / mass_model.W_scale;
save C:\isbe\dev\mass_model\weights\weights_sd weights_sd

% 2.
% Compute mean var weights
k_shape = length(mass_model.L_shape);
k_tex = length(mass_model.L_tex);
weights_var = [ k_shape / sum(mass_model.L_shape)...
                k_tex / sum(mass_model.L_tex)] * ...
                mass_model.L_scale;
save C:\isbe\dev\mass_model\weights\weights_var weights_var

% 3.
% Compute FD weights
%cd C:\isbe\dev\weights
load C:\isbe\dev\files\u_files.mat

r_idx = randsample(101, 50);
save C:\isbe\dev\mass_model\weights\shape_modes\r_idx r_idx;

shape_offsets = [-15 -10 -5 5 10 15];
scale_offsets = [-.015 -.0125 -.01 -.0075 -.005 -.0025...
    .0025 .005 .0075 .01 .0125 .015];
save C:\isbe\dev\mass_model\weights\offsets shape_offsets scale_offsets;

weights_fd_shape = zeros(k_shape, 1);
RMS_shapes = zeros(50, 6, k_shape);

for ii = 1:k_shape
    r_name = ['C:\isbe\dev\weights\shape_modes\s_mode', zerostr(ii, 2)];
    [RMS_shape] = ...
    calculate_weights2(u_files1, shape_offsets, [], 'shape_mode', ii,...
    'file_out', r_name, 'indices', r_idx,...
    'model', 'C:\isbe\dev\models\model_i500_50K.mat');
    weights_fd_shape(ii) = mean(mean(RMS_shape) ./ abs(shape_offsets));
    RMS_shapes(:,:,ii) = RMS_shape;
end
save C:\isbe\dev\mass_model\weights\RMS_shapes RMS_shapes;

r_name = 'C:\isbe\dev\weights\scale_mode';
[dummy, RMS_scale] = ...
calculate_weights2(u_files1, [], scale_offsets, 'shape_mode', ii,...
    'file_out', r_name, 'indices', r_idx,...
    'model', 'C:\isbe\dev\models\model_i500_50K.mat');            
weights_fd_scale = mean(mean(RMS_scale) ./ abs(scale_offsets));
weights_fd = [mean(weights_fd_shape) 1] / weights_fd_scale;
weights_fd_full = [weights_fd_shape' ones(1, k_tex) weights_fd_scale];
save C:\isbe\dev\mass_model\weights\weights_fd weights_fd*
%%
% 4.
% Compute optimised weights
cd C:\isbe\dev\weights
idx_opt = randsample(1:101, 50); save idx_opt idx_opt;
optimise_weights(mass_model, u_files1, [2e-3 5e-4], 'weights_opt', idx_opt);
%%
% if we have time: re-optimise using different random set/starting points
idx_opt1 = randsample(1:101, 25);
wso1 = [1.5e-3 0.1e-3];
[weights_opt1, output1] =...
    optimise_weights(mass_model, u_files1, wso1, 'weights_opt', idx_opt1);
save C:\isbe\dev\mass_model\weights\weights_opt1 idx_opt1 wso1 weights_opt1 output1;

idx_opt2 = randsample(1:101, 25);
wso2 = [1.5e-3 0.3e-3];
[weights_opt2, output2] =...
    optimise_weights(mass_model, u_files1, wso2, 'weights_opt', idx_opt2);
save C:\isbe\dev\mass_model\weights\weights_opt2 idx_opt2 wso2 weights_opt2 output2;
%
idx_opt3 = randsample(1:101, 25);
wso3 = [2.5e-3 0.1e-3];
[weights_opt3, output3] =...
    optimise_weights(mass_model, u_files1, wso3, 'weights_opt3', idx_opt3);
save C:\isbe\dev\mass_model\weights\weights_opt3 idx_opt3 wso3 weights_opt3 output3;

idx_opt4 = randsample(1:101, 25);
wso4 = [2.5e-3 0.3e-3];
[weights_opt4, output4] =...
    optimise_weights(mass_model, u_files1, wso4, 'weights_opt4', idx_opt4);
save C:\isbe\dev\mass_model\weights\weights_opt4 idx_opt4 wso4 weights_opt4 output4;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%
% Calculate leave one out errors
%
% Load weights
%cd C:\isbe\dev\weights
load C:\isbe\dev\mass_model\weights\weights_var
load C:\isbe\dev\mass_model\weights\weights_sd
load C:\isbe\dev\mass_model\weights\weights_fd;
% load weights_opt;
load C:\isbe\dev\mass_model\models\model_i500_50K.mat
load C:\isbe\dev\files\u_files.mat

[er_var.com er_var.ind] = ...
    model_errors_loo2(mass_model, u_files1, 'weights', weights_var, 'plot', 0);
save C:\isbe\dev\mass_model\weights\er_var er_var
[er_sd.com er_sd.ind] = ...
    model_errors_loo2(mass_model, u_files1, 'weights', weights_sd);
save C:\isbe\dev\mass_model\weights\er_sd er_sd
[er_fd.com er_fd.ind] = ...
    model_errors_loo2(mass_model, u_files1, 'weights', weights_fd);
save C:\isbe\dev\mass_model\weights\er_fd er_fd
%%
[er_opt.com er_opt.ind] = ...
    model_errors_loo2(mass_model, u_files1, 'weights', weights_opt);
save C:\isbe\dev\mass_model\weights\er_opt er_opt

[er_opt1.com er_opt1.ind] = ...
    model_errors_loo2(mass_model, u_files1, 'weights', weights_opt1);
save C:\isbe\dev\mass_model\weights\er_opt1 er_opt1

[er_opt2.com er_opt2.ind] = ...
    model_errors_loo2(mass_model, u_files1, 'weights', weights_opt2);
save C:\isbe\dev\mass_model\weights\er_opt2 er_opt2

[er_opt3.com er_opt3.ind] = ...
    model_errors_loo2(mass_model, u_files1, 'weights', weights_opt3);
save C:\isbe\dev\mass_model\weights\er_opt3 er_opt3

[er_opt4.com er_opt4.ind] = ...
    model_errors_loo2(mass_model, u_files1, 'weights', weights_opt4);
save C:\isbe\dev\mass_model\weights\er_opt4 er_opt4
%%
% Also try leave-one-out errors for full Finite Difference weights
[er_fd_full.com er_fd_full.ind] = ...
    model_errors_loo2(mass_model, u_files1, 'weights', weights_fd_full);
save C:\isbe\dev\mass_model\weights\er_fd_full er_fd_full

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute errors for map of weights
%%
load C:\isbe\dev\files\u_files.mat

shape_points = (0:0.5:5)*1e-3;
tex_points = (0:0.1:1)*1e-3;

n_s_pts = length(shape_points);
n_t_pts = length(tex_points);
error_map = zeros(n_s_pts, n_t_pts);

for ii = 2:n_s_pts
    for jj = 2:n_t_pts
        map_idx = randsample(1:101, 50);
        weights = [shape_points(ii) tex_points(jj)];
        [er.com er.ind] = ...
            model_errors_loo2(mass_model, u_files1, 'weights', weights, 'indices', map_idx);
        er_name = ['C:\isbe\dev\weights\map\er',...
            zerostr(ii, 2), zerostr(jj, 2)];
        error_map(ii, jj) = mean(er.com.weights);
        save(er_name, 'er', 'map_idx');
        clear weights er er_name;
    end
end

save error_map error_map shape_points tex_points
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Plot figures for report
cd C:\isbe\dev\weights
load weights_var
load weights_sd
load weights_fd;
load weights_opt;

%
% mean combined appearance vectors;
model_var = combine(mass_model, weights_var);
model_sd = combine(mass_model, weights_sd);
model_fd = combine(mass_model, weights_fd);
model_opt = combine(mass_model, weights_opt);

figure;
subplot(2,2,1)
bar(std(model_var.combined_data), 'b');
set(gca, 'XTick', [1 k_shape k_shape+k_tex]);
subplot(2,2,2)
bar(std(model_sd.combined_data), 'g');
set(gca, 'XTick', [1 k_shape k_shape+k_tex]);
subplot(2,2,3)
bar(std(model_fd.combined_data), 'r');
set(gca, 'XTick', [1 k_shape k_shape+k_tex]);
subplot(2,2,4)
bar(std(model_opt.combined_data), 'c');
set(gca, 'XTick', [1 k_shape k_shape+k_tex]);

%%
% loo errors figures
%
load er_var
load er_sd
load er_fd
load er_opt

%
% Boxplot of errors
figure;
boxplot([er_var.com.weights er_sd.com.weights er_fd.com.weights er_opt.com.weights]);
hold on;
plot([0.75, 1.25], [mean(er_var.com.weights), mean(er_var.com.weights)], 'g', 'LineWidth', 1.5);
plot([1.75, 2.25], [mean(er_sd.com.weights), mean(er_sd.com.weights)], 'g', 'LineWidth', 1.5);
plot([2.75, 3.25], [mean(er_fd.com.weights), mean(er_fd.com.weights)], 'g', 'LineWidth', 1.5);
plot([3.75, 4.25], [mean(er_opt.com.weights), mean(er_opt.com.weights)], 'g', 'LineWidth', 1.5);

plot([3.5, 4], [16 16], 'g', 'LineWidth', 1.5);
plot([3.5, 4], [15 15], 'r', 'LineWidth', 1.0);
% manually insert legend
%%
% Plot error map
load error_map

figure; hold on;
mesh(shape_points, tex_points, error_map);
plot3(weights_var(1), weights_var(2), mean(er_var.com.weights),...
    'bs', 'MarkerSize', 10);
plot3(weights_sd(1), weights_sd(2), mean(er_sd.com.weights),...
    'r^', 'MarkerSize', 10);
plot3(weights_fd(1), weights_fd(2), mean(er_fd.com.weights),...
    'mo', 'MarkerSize', 10);
% plot3(weights_opt(1), weights_opt(2), mean(er_opt.com.weights),...
%     'gd', 'MarkerSize', 10);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Figures for finite difference method
load RMS_shapes;

[N M k_shape] = size(RMS_shapes);

%%
%assume 0 displacemnet not computed +/- symmetrical shape displacements
RMS_modes = squeeze(mean(RMS_shapes));
RMS_modes = [RMS_modes(1:end/2,:); zeros(1, k_shape); RMS_modes(end/2+1:end,:)];

for ii = 1:k_shape; shape_leg{ii} = ['mode ', num2str(ii)]; end
legend(shape_leg);

figure; plot([shape_offsets(1:end/2), 0 , shape_offsets(end/2+1:end)], RMS_modes);
title('Mean texture displacement caused by peturbing shape modes by a given constant');
xlabel('Magnitude of shape mode displacement');
ylabel('Texture displacement');
set(gca, 'XTick', [shape_offsets(1:end/2), 0 , shape_offsets(end/2+1:end)]);

saveas(gcf, 'C:\isbe\dev\weights\figures\shape_displace1.fig');

%%
RMS_weights = zeros(N, k_shape);
for ii = 1:k_shape
    RMS_weights(:,ii) = mean(RMS_shapes(:,:,ii) ./...
        repmat(abs(shape_offsets), N, 1), 2);
end
sd = std(RMS_weights);
mm = mean(RMS_weights);

figure; hold on;

for ii = 1:k_shape
    
    plot([ii ii], [mm(ii)+2*sd(ii) mm(ii)-2*sd(ii)], 'bx', 'markersize', 5);
    plot([ii ii], [mm(ii)+2*sd(ii) mm(ii)-2*sd(ii)], 'b:', 'linewidth', 1.5);
    plot(ii, mm(ii), 'rx', 'MarkerSize', 10);
end

plot([0 k_shape+1], [mean(mm) mean(mm)], 'g');

title('Line plots of texture displacement per unit of shape parameter displacement for shape modes 1 to 29');
xlabel('Shape mode');
ylabel('Texture displacement');
set(gca, 'XTick', [1:k_shape]);
saveas(gcf, 'C:\isbe\dev\weights\figures\shape_displace2.fig');

%%
load C:\isbe\dev\weights\scale_mode.mat
RMS_scale = [RMS_scale(:,1:end/2), zeros(N, 1), RMS_scale(:,end/2+1:end)];
figure; plot([scale_offsets(1:end/2), 0 , scale_offsets(end/2+1:end)], RMS_scale');
title('Texture displacement for each mass caused by peturbing the scale mode by a given constant');
xlabel('Magnitude of scale mode displacement');
ylabel('Texture displacement');
set(gca, 'XTick', [scale_offsets(1:end/2), 0 , scale_offsets(end/2+1:end)]);
saveas(gcf, 'C:\isbe\dev\weights\figures\scale_displace.fig');

figure; plot([scale_offsets(1:end/2), 0 , scale_offsets(end/2+1:end)], mean(RMS_scale), '-rx');
title('Mean texture displacement caused by peturbing the scale mode by a given constant');
xlabel('Magnitude of scale mode displacement');
ylabel('Texture displacement');
set(gca, 'XTick', [scale_offsets(1:end/2), 0 , scale_offsets(end/2+1:end)]);
saveas(gcf, 'C:\isbe\dev\weights\figures\scale_displace_mean.fig');
%
figure; plot([scale_offsets(1:end/2), 0 , scale_offsets(end/2+1:end)], RMS_scale', ':'); hold on;
plot([scale_offsets(1:end/2), 0 , scale_offsets(end/2+1:end)], mean(RMS_scale), 'k', ...
    'LineWidth', 2.0);
title('Texture displacement for each mass caused by peturbing the scale mode by a given constant');
xlabel('Magnitude of scale mode displacement');
ylabel('Texture displacement');
set(gca, 'XTick', [scale_offsets(1:end/2), 0 , scale_offsets(end/2+1:end)]);
saveas(gcf, 'C:\isbe\dev\weights\figures\scale_displace_both.fig');
%%
%%%%
% correct error map you fool;

load C:\isbe\dev\files\u_files.mat

shape_points = [0.01 0.5:0.5:5]*1e-3;
tex_points = [0.01 0.1:0.1:1]*1e-3;

n_s_pts = length(shape_points);
n_t_pts = length(tex_points);
error_map = zeros(n_s_pts, n_t_pts);
%%
for ii = 11:n_s_pts
    for jj = 2:n_t_pts
        er_name = ['C:\isbe\dev\weights\map\er',...
            zerostr(ii, 2), zerostr(jj, 2)];
        load(er_name);
        new_idx = setdiff(1:101, map_idx); 
 
        weights = [shape_points(ii) tex_points(jj)];
        [er.com2 er.ind2] = ...
            model_errors_loo2(mass_model, u_files1, 'weights', weights, 'indices', new_idx);
        
        error_map(ii, jj) = mean([er.com.weights; er.com2.weights]);
        save(er_name, 'er', 'map_idx');
        clear weights er er_name;
    end
end

save error_map error_map shape_points tex_points
%%
weights = [shape_points(11) tex_points(1)];
[er.com2 er.ind2] = model_errors_loo2(mass_model, u_files1, 'weights', weights);
error_map(11, 1) = mean(er.com.weights);
save('C:\isbe\dev\weights\map\er1101', 'er');































