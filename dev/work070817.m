%%
idx = sort(randsample(49,6));
for ii = 1:6
    kk = idx(ii);
    load(['C:\isbe\dev\masses\left_out_models\model', zerostr(kk, 3)]);
    load(['C:\isbe\dev\masses\', m_files(kk).name]);
    leave_one_out_model(mass_model, mass);
    clear mass_model mass;
end
%%
idx = sort(randsample(179,6));
load C:\isbe\dev\mass_model
for ii = 1:6
    kk = idx(ii);
    load(['C:\isbe\dev\masses\', m_files(kk).name]);
    leave_one_out_model(mass_model, mass);
    clear mass;
end
%%
%load C:\isbe\dev\mass_model
mass_model2 = mass_model; clear mass_model 
idx = sort(randsample(179,3));
for ii = 1:3
    kk = idx(ii);
    load(['C:\isbe\dev\masses\left_out_models\model', zerostr(kk, 3)]);
    load(['C:\isbe\dev\masses\', m_files(kk+1).name]);
    leave_one_out_model(mass_model, mass);
    leave_one_out_model(mass_model2, mass);
    clear mass_model mass;
end
%%
for ii = 1:8
    kk = idx1(ii);
    load(['C:\isbe\dev\masses\', m_files(kk).name]);
    leave_one_out_model(mass_model, mass);
    clear mass;
end
%%
figure('WindowStyle', 'docked');
plot(zz1(1:40), zz3.indie(1:40), 'rx');
title('original shape vs first combined error');

figure('WindowStyle', 'docked');
plot(zz2(1:40), zz3.indie(1:40), 'rx');
title('original tex vs first combined error');

figure('WindowStyle', 'docked');
plot(er_shape_s(1:40)', zz3.combined(1:40), 'rx');
title('model shape vs final combined error');

figure('WindowStyle', 'docked');
plot(zz3.indie(1:40), zz3.combined(1:40), 'rx');
title('first combined error vs final combined error');
%%
figure('WindowStyle', 'docked');
plot(er_shape_s([1:40, 140:179])', er_combined, 'rx');
title('model tex vs final combined error');
hold on;

figure('WindowStyle', 'docked');
plot(er_tex_s([1:40, 140:179])', er_combined, 'rx');
title('model tex vs final combined error');
hold on;
%plot(zz3.tex(idx2(1:10))', zz3.combined(idx2(1:10)), 'gx');

figure('WindowStyle', 'docked');
plot(abs(er_scale_s([1:40, 140:179]))', er_combined, 'rx');
title('model scale vs final combined error');

figure('WindowStyle', 'docked');
plot(er_ss_s([1:40, 140:179])', er_combined, 'rx');
hold on;
%plot(er_ss_s(idx2(1:10))', zz3.combined(idx2(1:10)), 'gx');
title('scaled shape vs final combined error');

figure('WindowStyle', 'docked');
plot3(er_scale_s([1:40, 140:179])', er_tex_s([1:40, 140:179])', er_combined, 'rx');
title('scaled shape and tex vs final combined error');
%%
ll = size(mass_model.mean_shape_pl, 1);
for ii = 1:40
    load(['C:\isbe\dev\masses\', m_files(ii).name]);
    [r c] = size(mass.subtract_ROI);
    zz3.indie(ii) = zz3.indie(ii) * ll / (r*c);
    zz3.combined(ii) = zz3.combined(ii) * ll / (r*c);
end
%%
for ii = 1:40
    load(['C:\isbe\dev\masses\', m_files(ii).name]);
    [rows cols] = size(mass.subtract_ROI);
    bigness(ii) = rows*cols;
end
%%
clear;

load masses\m_files

idx = randsample(179,50);
r_files = m_files(idx);
calculate_weights('mass_model', r_files, 'masses\', 'weights\RMS_random');
%%
cd C:\isbe\dev;
[er_o] = model_errors2('mass_model', m_files, 'masses\', weights_o); save weights\er_o er_o;
[er_c] = model_errors2('mass_model', m_files, 'masses\', weights_c); save weights\er_c er_c;
[er_z] = model_errors2('mass_model', m_files, 'masses\', weights_z); save weights\er_z er_z;