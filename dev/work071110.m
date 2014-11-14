shape_points = [1e-4 (0.5:0.5:5)*1e-3];
tex_points = [ 1e-5 (0.1:0.1:1)*1e-3];

% n_s_pts = length(shape_points);
% n_t_pts = length(tex_points);
% error_map = zeros(n_s_pts, n_t_pts);

for ii = 1:n_s_pts
    for jj = 1
        map_idx = randsample(1:101, 50);
        weights = [shape_points(ii) tex_points(jj)];
        [er.com er.ind] = ...
            model_errors_loo2(mass_model, u_files1, 'weights', weights, 'indices', map_idx);
        er_name = ['C:\isbe\dev\weights\map\er',...
            num2str(ii), num2str(jj)];
        error_map(ii, jj) = mean(er.com.weights);
        save(er_name, 'er', 'map_idx');
        clear weights er er_name;
        pack;
    end
end

for ii = 1
    for jj = 2:n_t_pts
        map_idx = randsample(1:101, 50);
        weights = [shape_points(ii) tex_points(jj)];
        [er.com er.ind] = ...
            model_errors_loo2(mass_model, u_files1, 'weights', weights, 'indices', map_idx);
        er_name = ['C:\isbe\dev\weights\map\er',...
            num2str(ii), num2str(jj)];
        error_map(ii, jj) = mean(er.com.weights);
        save(er_name, 'er', 'map_idx');
        clear weights er er_name;
        pack;
    end
end

save error_map error_map shape_points tex_points
