function log_flow_error(estimated_flow, gt_flow, logfile, prefix)

if ~exist('logfile','var') || isempty(logfile), logfile = 1; end
if ~exist('prefix','var') || isempty(prefix), prefix = ''; end

gt_flow_dir = angle(gt_flow);
gt_flow_mag = abs(gt_flow);

est_flow_dir = angle(estimated_flow);
est_flow_mag = abs(estimated_flow);

diff_error = gt_flow - estimated_flow;
diff_error = abs(diff_error);
show_vec_stats(diff_error, 'Euclidean diff', logfile, '  ');

dir_error = angle(gt_flow_dir .* conj(est_flow_dir));
show_vec_stats(dir_error, 'Direction diff', logfile, '  ');

mag_error = gt_flow_mag - est_flow_mag;
show_vec_stats(mag_error, 'Magnitude diff', logfile, '  ');

fprintf(logfile, '\n');