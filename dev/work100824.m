matrix = 'V';
std = 'E';
weight = 'U';
%
args.num_samples = 5000;
args.feature_shape = 'rect';
args.bg_dir = 'C:\isbe\asymmetry_project\data\synthetic_backgrounds\smooth512\train\';
args.win_size = 3;
args.num_levels = 5;
args.do_max = 0;
args.bg_stem = [];
args.bg_ratio = 0;
args.fg_ratio = 0.025;
args.rotate = 0;

rand('twister', 10);
x = generate_dt_training_data(args);

[sOut, e] = g03aa(matrix, std, weight, x, int32(ones(1,size(x,2))), ones(1,size(x,2)), 0, int32(size(x,2)));
display(['Number of modes for win_size = 3: ' num2str(find(e(:,3)>0.99, 1)) ' out of ' num2str(size(e,1))]);
display(['Total variance for win_size = 3: ' num2str(sum(e(:,1)))]);
%
args.num_samples = 5000;
args.feature_shape = 'rect';
args.bg_dir = 'C:\isbe\asymmetry_project\data\synthetic_backgrounds\smooth512\train\';
args.win_size = 1;
args.num_levels = 5;
args.do_max = 0;
args.bg_stem = [];
args.bg_ratio = 0;
args.fg_ratio = 0.025;
args.rotate = 0;

rand('twister', 10);
x = generate_dt_training_data(args);

[sOut, e] = g03aa(matrix, std, weight, x, int32(ones(1,size(x,2))), ones(1,size(x,2)), 0, int32(size(x,2)));
display(['Number of modes for win_size = 1: ' num2str(find(e(:,3)>0.99, 1)) ' out of ' num2str(size(e,1))]);
display(['Total variance for win_size = 1: ' num2str(sum(e(:,1)))]);
%%
args.num_samples = 5000;
args.feature_shape = 'rect';
args.bg_dir = 'C:\isbe\asymmetry_project\data\synthetic_backgrounds\smooth512\train\';
args.win_size = 1;
args.num_levels = 5;
args.do_max = 0;
args.bg_stem = [];
args.bg_ratio = 0;
args.fg_ratio = 0.025;
args.rotate = 1;

rand('twister', 10);
x = generate_dt_training_data(args);

[sOut, e] = g03aa(matrix, std, weight, x, int32(ones(1,size(x,2))), ones(1,size(x,2)), 0, int32(size(x,2)));
display(['Number of modes for rotated data (win size 1): ' num2str(find(e(:,3)>0.99, 1)) ' out of ' num2str(size(e,1))]);
display(['Total variance for rotated data (win size 1): ' num2str(sum(e(:,1)))]);

args.num_samples = 5000;
args.feature_shape = 'rect';
args.bg_dir = 'C:\isbe\asymmetry_project\data\synthetic_backgrounds\smooth512\train\';
args.win_size = 3;
args.num_levels = 5;
args.do_max = 0;
args.bg_stem = [];
args.bg_ratio = 0;
args.fg_ratio = 0.025;
args.rotate = 1;

rand('twister', 10);
x = generate_dt_training_data(args);

[sOut, e] = g03aa(matrix, std, weight, x, int32(ones(1,size(x,2))), ones(1,size(x,2)), 0, int32(size(x,2)));
display(['Number of modes for rotated data (win size 3): ' num2str(find(e(:,3)>0.99, 1)) ' out of ' num2str(size(e,1))]);
display(['Total variance for rotated data (win size 3): ' num2str(sum(e(:,1)))]);
%%
args.num_samples = 5000;
args.feature_shape = 'clock';
args.bg_dir = 'C:\isbe\asymmetry_project\data\synthetic_backgrounds\smooth512\train\';
args.win_size = 1;
args.num_levels = 5;
args.do_max = 0;
args.bg_stem = [];
args.bg_ratio = 0;
args.fg_ratio = 0.025;
args.rotate = 0;

rand('twister', 10);
x = generate_dt_training_data(args);

[sOut, e] = g03aa(matrix, std, weight, x, int32(ones(1,size(x,2))), ones(1,size(x,2)), 0, int32(size(x,2)));
display(['Number of modes for clock data: ' num2str(find(e(:,3)>0.99, 1)) ' out of ' num2str(size(e,1))]);
display(['Total variance for clock data: ' num2str(sum(e(:,1)))]);
%
args.num_samples = 5000;
args.feature_shape = 'clock';
args.bg_dir = 'C:\isbe\asymmetry_project\data\synthetic_backgrounds\smooth512\train\';
args.win_size = 1;
args.num_levels = 5;
args.do_max = 0;
args.bg_stem = [];
args.bg_ratio = 0;
args.fg_ratio = 0.025;
args.rotate = 1;

rand('twister', 10);
x = generate_dt_training_data(args);

[sOut, e] = g03aa(matrix, std, weight, x, int32(ones(1,size(x,2))), ones(1,size(x,2)), 0, int32(size(x,2)));


display(['Number of modes for rotated clock data: ' num2str(find(e(:,3)>0.99, 1)) ' out of ' num2str(size(e,1))]);
display(['Total variance for rotated clock data: ' num2str(sum(e(:,1)))]);
%%
pca1 = u_load('Z:\asymmetry_project\data\misc\synthetic_lines_dt_clock_pca01.mat');
pca2 = u_load('Z:\asymmetry_project\data\misc\synthetic_lines_dt_clock_pca02.mat');
pca3 = u_load('Z:\asymmetry_project\data\misc\synthetic_lines_dt_clock_pca03.mat');
pca4 = u_load('Z:\asymmetry_project\data\misc\synthetic_lines_dt_clock_pca04.mat');
pca5 = u_load('Z:\asymmetry_project\data\misc\synthetic_lines_dt_clock_pca05.mat');

means = [pca1.mean' pca2.mean' pca3.mean' pca4.mean' pca5.mean'];
modes1 = [pca1.modes(:,1) pca2.modes(:,1) pca3.modes(:,1) pca4.modes(:,1) pca5.modes(:,1)];

mean_errs = std(means,0, 2);
modes_errs_modes = std(modes1,0, 2);

save('Z:\asymmetry_project\data\misc\synthetic_lines_dt_clock_pca_errs_new.mat', '*errs');
%%