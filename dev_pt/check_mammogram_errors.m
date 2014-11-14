% clc;
% clear;

mampath = [asymmetryroot('shared'),'data\synthetic_lines\real512\'];

% default filename format (overwritten where need be)
format_str = '%02d_test_ext_class.mat';

respath = [mampath,'predictions\g2d\analytic\']; format_str = 'orientations/%03d_test_ori.mat';

orientations_vec = [];
for i = 1:100
    est_ori = load_uint8([respath,sprintf(format_str,i)]); % estimated
	load([mampath,'masks\',sprintf('mask%03d.mat',i)]); % ground truth
    orientations_vec = [orientations_vec; est_ori(mask)];
end

% debug: compare with errors estimated by synth_mammo_tests.m
if ~exist('gt_orientations','var')
	load([mampath,'orientations\all_gt_orientations.mat']);
end
errors_vec = ori_error(orientations_vec,gt_orientations);
median(abs(errors_vec) * 180/pi)

if ~exist('prediction_errs','var')
	load([respath,'errors/orientation_errors.mat']);
end
median(abs(prediction_errs) * 180/pi)

