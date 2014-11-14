f_com = figure('WindowStyle', 'normal', 'Position', [700,450, 600, 600]);
f_tex = figure('WindowStyle', 'normal', 'Position', [100,450, 600, 600]);
f_shape = figure('WindowStyle', 'normal', 'Position', [400,50, 400, 400]);

movie(f_shape, shape_frames, -3);
movie(f_tex, tex_frames, -3);
movie(f_com, com_frames, -3);

%%

load C:\isbe\dev\files\bg_files.mat
%%

gd = [2 6 7 11 12 16 19 23 24 27 28 32 33 34 35 36 39 40 43 45];
for ii = 1:length(gd)
    load(['C:\isbe\dev\masses\', bg_files(gd(ii)).name]);
    figure('name', ['Mass ', num2str(gd(ii))]);
    cmax = max(mass.mass_ROI(:));
    cmin = min(mass.mass_ROI(:));
    subplot(2,2,1); imagesc(mass.mass_ROI); axis image; colormap(gray(256));
    subplot(2,2,2); imagesc(double(mass.mass_ROI) - mass.subtract_ROI);
    axis image; colormap(gray(256)); caxis([cmin cmax]);
    subplot(2,2,3:4); imagesc(mass.subtract_ROI);
    axis image; colormap(gray(256)); hold on;
    %plot(mass.mass_outline(:,1), mass.mass_outline(:,2), 'g:');
    clear mass;
end
%%
clear
load C:\isbe\dev\files\u_files.mat
load C:\isbe\dev\models\model_w500_50K.mat
%%
model_errors_loo2(mass_model, u_files1, 'indices', [21:40]);
%%
load(['C:\isbe\dev\masses\', bg_files(28).name]);
[i_error c_error new_mass] = leave_one_out_model(mass_model, mass, bg_good(28));
%%
load(['C:\isbe\dev\masses\', bg_files(39).name]);
[i_error c_error new_mass] = leave_one_out_model(mass_model, mass, bg_good(39));
%%
load(['C:\isbe\dev\masses\', bg_files(43).name]);
[i_error c_error new_mass] = leave_one_out_model(mass_model, mass, bg_good(43));