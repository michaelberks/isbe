%Work 20090601 - June 2009!!!!!!!!!!!!!!!!!!!!!!!!!!
normal_list = dir('C:\isbe\dev\background\images\normal1024\*.bmp');
[mass_model model_id] = u_load('C:\isbe\dev\mass_model\models\model_w500_50K.mat');
%
mass_idx = randperm(length(model_id.mass_files));
mass_dir = 'C:\isbe\dev\masses1024x1024';

for ii = 1:length(normal_list);
    target_region = double(imread(['C:\isbe\dev\background\images\normal1024\', normal_list(ii).name]));
    for ra = 0:30:330       
        rotation_angle = pi*ra/180;
        mass = mb_synthesise_mass(target_region, [512 512], mass_model, mass_dir, model_id.mass_files, mass_idx(ii), rotation_angle);
        mass_name = ['C:\isbe\dev\synthetic_masses\n', zerostr(ii,3), '_m', zerostr(mass_idx(ii),3), '_a', zerostr(ra,3), '.bmp'];
        imwrite(uint8(mass.mass_ROI), mass_name);
    end
end
%%
norm_idx = randperm(89);
for ii = 1:30
    copyfile(['C:\isbe\dev\synthetic_masses\n', zerostr(norm_idx(ii),3), '*'], 'C:\isbe\dev\observer_study\masses\')
end
%%
mass_list = dir('C:\isbe\dev\observer_study\masses\*.mat');
for ii = 1:length(mass_list);
    load(['C:\isbe\dev\observer_study\masses\', mass_list(ii).name]);
    imwrite(uint8(mass_ROI), ['C:\isbe\dev\observer_study\masses\', mass_list(ii).name(1:end-3), 'bmp']);
end

%%
normal_list = dir('C:\isbe\dev\background\images\normal1024\*.bmp');
[mass_model model_id] = u_load('C:\isbe\dev\mass_model\models\model_w500_50K.mat');

for ii = 1:30;
    mass_idx = template_mass_idx(ii);
    
    sample_args.TargetRegion = double(imread(['C:\isbe\dev\background\images\normal1024\', normal_list(norm_idx(ii)).name]));
    sample_args.TargetCentre = [512 512];
    sample_args.ConditionMode = mass_idx;
    sample_args.MassModel = mass_model;
    sample_args.NumModes = 5;
    rotation_angle = pi*angles(ii)/180;
    sample_args.Rotation = mass_model.rotations(:,:,mass_idx)*...
        [cos(rotation_angle) sin(rotation_angle); -sin(rotation_angle) cos(rotation_angle)];
    sample_args.Origin = mass_model.origins(mass_idx);
    
    [mass] = sample_new_mass_in_region(sample_args);
    mass.template_idx = mass_idx;
    mass_name = ['C:\isbe\dev\synthetic_masses\n', zerostr(norm_idx(ii),3), '_m', zerostr(mass_idx,3), '_a', zerostr(angles(ii),3), 'no_bg.bmp'];
    imwrite(uint8(mass.mass_ROI), mass_name);

end
%%
normal_list = dir('C:\isbe\dev\background\images\normal1024\*.bmp');
target_region = (double(imread(['C:\isbe\dev\background\images\normal1024\', normal_list(24).name])));
% tree = imread('C:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\tree\tree_gray.bmp');
% tree = tree(1:1024, 1:1024);
%build dual-tree
target_dt = dtwavexfm2(target_region, 6);

%convert to the full tree
[target_ilp target_icp] = mb_dual_tree_transform(target_dt);
clear target_dt target_ilp
max_icp = cell(4,1);
max_icp_full = zeros(1024,1024,4);
for lev = 1:4
    max_icp{lev} = max(target_icp{lev+1}, [], 3);
    max_icp_full(:,:,lev) = kron(max_icp{lev}, ones(2^(lev+1)));
end
clear target_icp
%
figure; image(complex2rgb(max_icp_full(:,:,4), [0 pi])); axis image; hold on;
figure; image(complex2rgb(max_icp_full(:,:,3), [0 pi])); axis image; hold on;
figure; image(complex2rgb(max_icp_full(:,:,2), [0 pi])); axis image; hold on;
figure; image(complex2rgb(max_icp_full(:,:,1), [0 pi])); axis image; hold on;
%
figure; weighted_complex_rose(max_icp_full(:,:,1), 100);
figure; weighted_complex_rose(max_icp_full(:,:,2), 100);
figure; weighted_complex_rose(max_icp_full(:,:,3), 100);
figure; weighted_complex_rose(max_icp_full(:,:,4), 100);
%
max_icp_full_all = max(max_icp_full, [], 3);
figure; image(complex2rgb(max_icp_full_all, [0 pi])); axis image; hold on;
figure; weighted_complex_rose(max_icp_full_all, 360);
[dummy dummy bins] = weighted_complex_rose(max_icp_full_all, 360);
%%
[orientation_map] = mb_dt_orientation_map('Image', target_region);
