function student_symposium_script

% Masses I like

ok_masses = [8 11 14 15 17 21 23 29 32 60 68 72 73 74 77 84 85 88 89 91 98 100];
%%
for ii = 1:length(ok_masses)
    load(['C:\isbe\dev\masses\', u_files1(ok_masses(ii)).name]);
    figure; imagesc(mass.subtract_ROI); colormap(gray(256)); axis image;
    clear mass;
end

load C:\isbe\dev\models\model_o500_50K
%%
for idx = [8 19]
    load(['C:\isbe\dev\masses\', u_files1(ok_masses(idx)).name]);
    leave_one_out_model(mass_model, mass, ok_masses(idx));
    clear mass
end
%%
for ii = 1:length(ok_masses)
    load(['C:\isbe\dev\masses\', u_files1(ok_masses(ii)).name]);
    write_im_from_colormap(mass.subtract_ROI, ['C:\isbe\ss08\ok_mass', zerostr(ii,3), '.bmp']);
    figure; imagesc(mass.subtract_ROI); colormap(gray(256)); axis image;
    clear mass;
end