cd C:\isbe\dev\models
%%
load C:\isbe\dev\models\origin_model.mat
generate_loo_models('mass_model', mass_model, 'model_path', 'loo_u1\');
clear mass_model;
load C:\isbe\dev\models\small_model.mat
generate_loo_models('mass_model', mass_model, 'model_path', 'loo_small\');
clear mass_model;
load C:\isbe\dev\models\large_model.mat
generate_loo_models('mass_model', mass_model, 'model_path', 'loo_large\');
clear mass_model;
load C:\isbe\dev\models\r51_model.mat
generate_loo_models('mass_model', mass_model, 'model_path', 'loo_r51\');
clear mass_model;
load C:\isbe\dev\models\r50_model.mat
generate_loo_models('mass_model', mass_model, 'model_path', 'loo_r50\');
clear mass_model;
%%
load C:\isbe\dev\files\u_files.mat
cd C:\isbe\dev\weights\shape_modes
r_idx = randsample(101, 50);
save r_idx r_idx;

for ii = 2:31
    r_name = ['C:\isbe\dev\weights\s_mode', zerostr(ii, 2)];
    calculate_weights2(u_files1, [-15 -10 -5 5 10 15], [], 'shape_mode', ii,...
    'file_out', r_name, 'indices', r_idx,...
    'model', 'C:\isbe\dev\models\origin_model.mat');
end
calculate_weights2(u_files1, [], [-.015 -.01 -.005 .005 .01 .015], 'shape_mode', ii,...
    'file_out', r_name, 'indices', r_idx,...
    'model', 'C:\isbe\dev\models\origin_model.mat');