load C:\isbe\dev\files\u_files.mat
idx = randsample(101, 20);

for ii = 1:20
    jj = idx(ii);
    load(['C:\isbe\dev\masses\', u_files1(jj).name]);
    %mass_true = mass.subtract_ROI > 0;
    figure; hold on;
    %imagesc(mass_true);
    imagesc(mass.subtract_ROI); colormap(gray(256)); axis image;
    plot(mass.dilate_outline(:,1), mass.dilate_outline(:,2), 'y');
    clear mass mass_true;
end

%%

idx = randsample(101, 20);

for ii = 1:20
    jj = idx(ii);
    load(['C:\isbe\dev\masses\', u_files1(jj).name]);
    figure; hold on;
    imagesc(mass.subtract_ROI); colormap(gray(256)); axis image;
    plot(mass.dilate_outline(:,1), mass.dilate_outline(:,2), 'y');
    plot(mass.mass_outline(:,1) + mass.mass_centroid(1), mass.mass_outline(:,2) + mass.mass_centroid(2), 'b');
    clear mass mass_true;
    
    load(['C:\isbe\dev\annotations\', u_files1(jj).name]);
    plot(mass.mass_outline(:,1) + mass.mass_centroid(1), mass.mass_outline(:,2) + mass.mass_centroid(2), 'r:');
    plot(mass.annotated_outline(:,1), mass.annotated_outline(:,2), 'g:');
end
%%
load C:\isbe\dev\files\large_files.mat
load C:\isbe\dev\files\small_files.mat
load C:\isbe\dev\files\r_files50.mat
load C:\isbe\dev\files\r_files51.mat

generate_mass_AM(large_files, 'C:\isbe\dev\models\large_model', 500, 50000);
generate_mass_AM(small_files, 'C:\isbe\dev\models\small_model', 500, 50000);
generate_mass_AM(r_files50, 'C:\isbe\dev\models\r50_model', 500, 50000);
generate_mass_AM(r_files51, 'C:\isbe\dev\models\r51_model', 500, 50000);
clear

%%
load C:\isbe\dev\files\large_files.mat
cd C:\isbe\dev\models\loo_large\

for ii = 1:51
    model_name = ['model', zerostr(ii, 3)];
    load(model_name);
    loo_files = large_files;
    loo_files(ii) = [];
    change_model_shape(mass_model, model_name, loo_files, 500)
    clear mass_model
end
clear
%%
load C:\isbe\dev\files\small_files.mat
cd C:\isbe\dev\models\loo_small\

for ii = 1:50
    model_name = ['model', zerostr(ii, 3)];
    load(model_name);
    loo_files = small_files;
    loo_files(ii) = [];
    change_model_shape(mass_model, model_name, loo_files, 500)
    clear mass_model
end
%%
load C:\isbe\dev\files\r_files50.mat
cd C:\isbe\dev\models\loo_r50\

for ii = 1:50
    model_name = ['model', zerostr(ii, 3)];
    load(model_name);
    loo_files = r_files50;
    loo_files(ii) = [];
    change_model_shape(mass_model, model_name, loo_files, 500)
    clear mass_model
end
%%
load C:\isbe\dev\files\u_files.mat
cd C:\isbe\dev\models\loo_u1\

for ii = 1:101
    model_name = ['model', zerostr(ii, 3)];
    load(model_name);
    loo_files = u_files1;
    loo_files(ii) = [];
    change_model_shape(mass_model, model_name, loo_files, 500)
    clear mass_model
end
%
load C:\isbe\dev\files\r_files51.mat
cd C:\isbe\dev\models\loo_r51\

for ii = 1:41
    model_name = ['model', zerostr(ii, 3)];
    load(model_name);
    loo_files = r_files51;
    loo_files(ii) = [];
    change_model_shape(mass_model, model_name, loo_files, 500)
    clear mass_model
end