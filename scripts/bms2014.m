fig_dir = 'K:\isbe\conferences_and_symposia\BMS2014\figs\';

nailfold = u_load('C:\isbe\nailfold\data\rsa_study\master_set\images\10598c.mat');
vessel_prob = u_load('C:\isbe\nailfold\data\rsa_study\master_set\predictions\detection\rf_classification\296655\10598c_pred.mat');
vessel_ori = u_load('C:\isbe\nailfold\data\rsa_study\master_set\predictions\orientation\rf_regression\296621\10598c_pred.mat');
vessel_centre = u_load('C:\isbe\nailfold\data\rsa_study\master_set\vessel_centres\full_centres\10598c_vc.mat');
load('C:\isbe\nailfold\data\rsa_study\master_set\apex_maps\set12g_half_296655\10598c_pred.mat');
load('C:\isbe\nailfold\data\rsa_study\master_set\apex_maps\set12g_half_296655\miccai_maxima\rescores\10598c_candidates.mat');
load('C:\isbe\nailfold\data\rsa_study\master_set\apex_maps\set12g_half_296655\miccai_maxima\selected_apexes\10598c_sel.mat');
load('C:\isbe\nailfold\data\rsa_study\master_set\apex_clusters_merged\10598c_apex_clusters.mat');

initial_thresh = 0.3;
apex_class_thresh = 0.5;

cols = 267:495;
rows = 160:280;
ncols = length(cols);
nrows = length(rows);

vessel_centre.x = vessel_centre.x - cols(1) + 1;
vessel_centre.y = vessel_centre.y - rows(1) + 1;
candidate_xy(:,1) = candidate_xy(:,1) - cols(1) + 1;
candidate_xy(:,2) = candidate_xy(:,2) - rows(1) + 1;

nailfold = nailfold(rows, cols);
vessel_prob = vessel_prob(rows, cols);
vessel_ori = vessel_ori(rows, cols);
apex_offset_map = apex_offset_map(rows,cols);
include_pts = apex_class_pred > apex_class_thresh;

valid_candidates = candidate_rescores > initial_thresh;

xx = 1:ncols;
yy = (1:nrows)';
grid_x = repmat(xx, nrows, 1);
grid_y = repmat(yy, 1, ncols);
   
%Compute weighted kernel estimates of the spatial distribution of
%candidates over this grid
[location_distribution] = build_2d_kernel_distribution(...
   candidate_xy(valid_candidates,:),...
   [grid_x(:) grid_y(:)],...
   candidate_rescores(valid_candidates,:), 0);
   
location_distribution.D_f = reshape(location_distribution.D_f, size(grid_x));
[~, y_max] = max(location_distribution.D_f);
candidate_polyfit = interp1(xx, yy(y_max), candidate_xy(:,1), 'linear');

%%
figure; imgray(vessel_prob);
axis off;
exportfig_im([fig_dir '01_vessel_prediction.png'], [nrows ncols]);
plot(vessel_centre.x, vessel_centre.y, 'g.', 'markersize', 8);
exportfig_im([fig_dir '01_vessel_prediction_vc.png'], [nrows ncols]);
plot(vessel_centre.x(include_pts), vessel_centre.y(include_pts), 'r.', 'markersize', 10);
exportfig_im([fig_dir '01_vessel_prediction_vc2.png'], [nrows ncols]);

%
figure; imgray(complex2rgb(vessel_ori));
axis off;
exportfig_im([fig_dir '02_orientation_prediction.png'], [nrows ncols]);
%
figure; imgray(apex_offset_map); colormap(hot(256));
axis off;
exportfig_im([fig_dir '03_apex_heat_map.png'], [nrows ncols]);
%
figure; imgray(location_distribution.D_f);
axis off;
plot(1:length(cols), y_max, 'g--', 'linewidth', 2);
for i_c = 1:length(candidate_rescores)
    if candidate_xy(i_c,1) > 1 && candidate_xy(i_c,1) < ncols &&...
            candidate_xy(i_c,2) > 1 && candidate_xy(i_c,2) < nrows
    plot(candidate_xy(i_c,1), candidate_xy(i_c,2), 'ro', 'markersize', 10*candidate_rescores(i_c));
    plot(...
        [candidate_xy(i_c,1) candidate_xy(i_c,1)],...
        [candidate_xy(i_c,2) candidate_polyfit(i_c)], 'r');
    end
end
exportfig_im([fig_dir '04_apex_location_map.png'], [nrows ncols]);