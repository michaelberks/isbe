for ii = 1:179
    load(['C:\isbe\dev\annotations\', anno_list(ii).name]);
    
    ROI = uint8(double(mass.mass_ROI) - double(mass.mass_spline));
    clear mass;
    
    load(['C:\isbe\dev\sjc_masses\', m_files(ii).name]);
    mass.subtract_ROI = ROI;
    save(['C:\isbe\dev\sjc_masses\', m_files(ii).name], 'mass');
    clear mass;
    
end

%%
r_idx = randsample(179, 20);
for ii = 1:20
    figure;
    kk = r_idx(ii);
    load(['C:\isbe\dev\sjc_masses\', m_files(ii).name]);
    subplot(1,2,1); imagesc(mass.subtract_ROI); axis image; colormap(gray(256));
    clear mass;
    
    load(['C:\isbe\dev\masses\', m_files(ii).name]);
    subplot(1,2,2); imagesc(mass.subtract_ROI); axis image; colormap(gray(256));
    clear mass;
end
    
%%
for ii = 1:179
    
    load(['C:\isbe\dev\sjc_masses\', m_files(ii).name]);
    mass.subtract_ROI = double(mass.subtract_ROI);
    save(['C:\isbe\dev\sjc_masses\', m_files(ii).name], 'mass');
end

%%
load er_var
load er_sd
load er_fd
load er_opt
load er_opt1
load er_opt2
load er_opt3
load er_opt4
load er_fd_full

%
% Boxplot of errors
figure;
boxplot([er_var.com.weights er_sd.com.weights er_fd.com.weights er_opt.com.weights er_sjc.com.weights er_fd_full.com.weights]);
hold on;
plot([0.75, 1.25], [mean(er_var.com.weights), mean(er_var.com.weights)], 'g', 'LineWidth', 1.5);
plot([1.75, 2.25], [mean(er_sd.com.weights), mean(er_sd.com.weights)], 'g', 'LineWidth', 1.5);
plot([2.75, 3.25], [mean(er_fd.com.weights), mean(er_fd.com.weights)], 'g', 'LineWidth', 1.5);
plot([3.75, 4.25], [mean(er_opt.com.weights), mean(er_opt.com.weights)], 'g', 'LineWidth', 1.5);
plot([4.75, 5.25], [mean(er_sjc.com.weights), mean(er_sjc.com.weights)], 'g', 'LineWidth', 1.5);  
plot([5.75, 6.25], [mean(er_fd_full.com.weights), mean(er_fd_full.com.weights)], 'g', 'LineWidth', 1.5);

mean(er_var.com.weights)
mean(er_sd.com.weights)
mean(er_fd.com.weights)
mean(er_opt.com.weights)
mean(er_opt1.com.weights)
mean(er_opt2.com.weights)
mean(er_opt3.com.weights)
mean(er_opt4.com.weights)
mean(er_sjc.com.weights)
mean(er_fd_full.com.weights)

mean(er_var.com.total)
mean(er_sd.com.total)
mean(er_fd.com.total)
mean(er_opt.com.total)
mean(er_opt1.com.total)
mean(er_opt2.com.total)
mean(er_opt3.com.total)
mean(er_opt4.com.total)
mean(er_sjc.com.total)
mean(er_fd_full.com.total)

%%
dirName = 'C:\isbe\dev\mdl\isbi';
    
amb_automodelbuild (reshape_MDL(shapes_aligned), 'circle', 'nIterations', Inf,...
    'saveFrequency', 50, 'Quiet', 1,...
    'optimisePose', 1, 'optimiseOrigin', 1, 'nExamplesToOptimise', 101,...
    'initialAlign', 0, 'initialOriginOptimisation', 0,...
    'nSamplePoints', 500, 'totalMaxFevals', 1e5, 'saveDirName', dirName);
%%
load C:\isbe\dev\weights\error_map
load C:\isbe\dev\weights\er_sd
load C:\isbe\dev\weights\er_fd
load C:\isbe\dev\weights\er_var
load C:\isbe\dev\weights\er_opt3


fine_shape_points = (0:0.25:5)*1e-3;
fine_tex_points = (0:0.05:1)*1e-3;

shape_points = (0:0.5:5)*1e-3;
tex_points = (0:0.1:1)*1e-3;

[x_pts y_pts] = meshgrid(shape_points, tex_points);
[fx_pts fy_pts] = meshgrid(fine_shape_points, fine_tex_points);
fine_map = interp2(x_pts, y_pts, error_map, fx_pts, fy_pts);

figure; hold on;
mesh(fine_shape_points, fine_tex_points, fine_map);
xlabel('W_{shape}');
ylabel('W_{tex}');
hidden on;

load C:\isbe\dev\weights\weights_sd
load C:\isbe\dev\weights\weights_fd
load C:\isbe\dev\weights\weights_var
load C:\isbe\dev\weights\weights_opt3
%%
var_pt = interp2(x_pts, y_pts, error_map, weights_var(1), weights_var(2));
sd_pt = interp2(x_pts, y_pts, error_map, weights_sd(1), weights_sd(2));
fd_pt = interp2(x_pts, y_pts, error_map, weights_fd(1), weights_fd(2));
opt_pt = interp2(x_pts, y_pts, error_map, weights_opt3(1), weights_opt3(2));

plot3(weights_var(1), weights_var(2), var_pt,...
    'bs', 'MarkerSize', 10);
plot3(weights_sd(1), weights_sd(2), sd_pt,...
    'r^', 'MarkerSize', 10);
plot3(weights_fd(1), weights_fd(2), fd_pt,...
    'mo', 'MarkerSize', 10);
plot3(weights_opt3(1), weights_opt3(2), opt_pt,...
     'gd', 'MarkerSize', 10);
 
%%
% Display our model modes