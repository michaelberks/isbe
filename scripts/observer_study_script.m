%Script for putting together the dataset needed for the observer study

%Step1 is to convert our original masses into structures containing
%1024x1024 regions
mass_list = dir('C:\isbe\dev\masses\*.mat');
mb_change_mass_region_size('MassList', mass_list, 'SavePath', 'C:\isbe\dev\masses1024x1024\');

%Now we need to select 30 random masses from the real masses not used in
%building the model
load('C:\isbe\dev\files\u_files.mat');
unused_masses_idx = setdiff(1:179, idx_u1);

%Now take a random ordering of these masses
unused_masses_idx = unused_masses_idx(randperm(length(unused_masses_idx)));
%
%Now take the first 30 of these that have full 1024x1024 regions
ii = 1;
jj = 1;
real_mass_idx = zeros(30,1);
while ii <= 30
    load(['C:\isbe\dev\masses1024x1024\', mass_list(unused_masses_idx(jj)).name])
    
    if all(size(mass.mass_ROI) == [1024 1024])
        mass_ROI = mass.mass_ROI;
        figure; imagesc(mass_ROI); colormap(gray(256)); axis image;
        real_mass_idx(ii) = unused_masses_idx(jj);
        save(['C:\isbe\dev\observer_study\masses\real_mass', zerostr(ii,1)], 'mass_ROI');
        ii = ii+1;
    end
    jj = jj+1;
end
save('C:\isbe\dev\observer_study\real_mass_idx', 'real_mass_idx');

%%
%load('C:\isbe\dev\files\u_files.mat');
%Now create the synthetic masses...
for ii = 1:11;
    mass_real = u_load(['C:\isbe\dev\masses1024x1024\', u_files1(ii).name]);
    args.TargetRegion = mass_real.background_ROI;
    args.TargetCentre = mean(mass_real.mass_outline);
    args.ConditionMode = ii;
    args.Rotation = mass_model.rotations(:,:,ii);
    args.Origin = mass_model.origins(ii);
    load C:\isbe\dev\mass_model\models\model_w500_50K.mat
    args.MassModel = mass_model;
    args.Plot = 0;
    args.NumModes = 32;
    [mass] = sample_new_mass_in_region(args);
    figure;
    subplot(1,2,1); imagesc(mass_real.mass_ROI - mass_real.background_ROI); axis image; colormap(gray(256));
    subplot(1,2,2); imagesc(mass.subtract_ROI); axis image; colormap(gray(256));
    %save(['C:\isbe\dev\background\new_masses\real_copy5\mass', zerostr(ii,3)], 'mass');
end
%%
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
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
syn_a_list = dir('C:\isbe\dev\observer_study\masses\bmp\syn_A\*.bmp');
syn_b_list = dir('C:\isbe\dev\observer_study\masses\bmp\syn_B\*.bmp');

for ii = 1:30
    mass_ROI = double(imread(['C:\isbe\dev\observer_study\masses\bmp\syn_A\', syn_a_list(ii).name]));
    save(['C:\isbe\dev\observer_study\masses\syn_a_mass', zerostr(ii,2)], 'mass_ROI');
    
    mass_ROI = double(imread(['C:\isbe\dev\observer_study\masses\bmp\syn_B\', syn_b_list(ii).name]));
    save(['C:\isbe\dev\observer_study\masses\syn_b_mass', zerostr(ii,2)], 'mass_ROI');
end
%%
syn_c_list = dir('C:\isbe\dev\observer_study\masses\bmp\syn_C\*.bmp');

for ii = 1:30
    mass_ROI = double(imread(['C:\isbe\dev\observer_study\masses\bmp\syn_C\', syn_c_list(ii).name]));
    save(['C:\isbe\dev\observer_study\masses\syn_c_mass', zerostr(ii,2)], 'mass_ROI');
end
%%
%Work out which normal regions and mass templates we used we used
norm_idx = zeros(30,1);
mass_idx = zeros(30,1);
syn_list = dir('C:\isbe\dev\observer_study\masses\bmp\syn_A\*.bmp');
for ii = 1:30
    norm_idx(ii) = str2num(syn_list(ii).name(2:4));
    mass_idx(ii) = str2num(syn_list(ii).name(7:9));
end
%
normal_list = dir('C:\isbe\dev\background\images\normal1024\*.bmp');
[mass_model model_id] = u_load('C:\isbe\dev\mass_model\models\model_w500_50K.mat');
mass_dir = 'C:\isbe\dev\masses1024x1024';
%%
for ii = 1:30
    target_region = double(imread(['C:\isbe\dev\background\images\normal1024\', normal_list(norm_idx(ii)).name]));
    [mass ra] = ...
        mb_synthesise_mass(target_region, [512 512], mass_model,...
        mass_dir, model_id.mass_files, mass_idx(ii), []);
    ra = round(180*ra/pi);
    mass_name = ['C:\isbe\dev\observer_study\masses\bmp\syn_C\n', ...
        zerostr(norm_idx(ii),3), '_m', zerostr(mass_idx(ii),3), '_a', zerostr(ra,3), '_c', '.bmp'];
    imwrite(uint8(mass.mass_ROI), mass_name);
end
%%
syn_a_list = dir('C:\isbe\dev\observer_study\masses\bmp\syn_A\*.bmp');
for ii = 1:30
    load(['C:\isbe\dev\observer_study\masses\syn_c_mass', zerostr(ii,2)], 'mass_ROI');
    imwrite(uint8(mass_ROI), ['C:\isbe\dev\observer_study\masses\bmp\syn_B\', syn_a_list(ii).name]);
end