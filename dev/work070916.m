%%
load C:\isbe\dev\files\u_files.mat
for ii = 21:101
    file_out = ['C:\isbe\dev\models\loo\model', zerostr(ii, 3)];
    one_out_files = u_files1;
    one_out_files(ii) = [];

    mass_model = generate_mass_AM(...
        one_out_files, file_out, 500, 'C:\isbe\dev\masses\', 20000);
end
%%
load C:\isbe\dev\files\large_files.mat
for ii = 1:51
    file_out = ['C:\isbe\dev\models\loo_large\model', zerostr(ii, 3)];
    one_out_files = large_files;
    one_out_files(ii) = [];

    mass_model = generate_mass_AM(...
        one_out_files, file_out, 500, 'C:\isbe\dev\masses\', 20000);
end

load C:\isbe\dev\files\small_files.mat
for ii = 1:50
    file_out = ['C:\isbe\dev\models\loo_small\model', zerostr(ii, 3)];
    one_out_files = small_files;
    one_out_files(ii) = [];

    mass_model = generate_mass_AM(...
        one_out_files, file_out, 500, 'C:\isbe\dev\masses\', 20000);
end
%%
er_small = model_errors_loo(small_files, weights_u, 'C:\isbe\dev\models\loo_small\model');
save C:\isbe\dev\scale\er_small er_small;
er_large = model_errors_loo(large_files, weights_u, 'C:\isbe\dev\models\loo_large\model');
save C:\isbe\dev\scale\er_large er_large;

for ii = 1:51
    file_out = ['C:\isbe\dev\models\loo_r51\model', zerostr(ii, 3)];
    one_out_files = r_files51;
    one_out_files(ii) = [];

    mass_model = generate_mass_AM(...
        one_out_files, file_out, 500, 'C:\isbe\dev\masses\', 20000);
end
%%
load C:\isbe\dev\scale\weights_u.mat
load C:\isbe\dev\files\small_files.mat
load C:\isbe\dev\files\large_files.mat
load C:\isbe\dev\files\r_files50.mat
load C:\isbe\dev\files\r_files51.mat

%er_u = model_errors3(u_files1, weights_u);
%save C:\isbe\dev\scale\er_u er_u;

er_small = model_errors_loo(small_files, weights_u, 'C:\isbe\dev\models\loo_small');
save C:\isbe\dev\scale\er_small er_small;
er_large = model_errors_loo(large_files, weights_u, 'C:\isbe\dev\models\loo_small');
save C:\isbe\dev\scale\er_large er_large;
er_r50 = model_errors_loo(r_files50, weights_u, 'C:\isbe\dev\models\loo_r50');
save C:\isbe\dev\scale\er_r50 er_r50;
er_51 = model_errors_loo(r_files51, weights_u, 'C:\isbe\dev\models\loo_r51');
save C:\isbe\dev\scale\er_51 er_51;