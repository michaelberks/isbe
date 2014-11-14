% check Staal's Az

clear; clc;

aroot = asymmetryroot;
img_path = [aroot, 'data\retinograms\drive\test\images\'];
staal_path = [aroot, 'data\retinograms\drive\test\output\Staal_etal\'];
msk_path = [aroot, 'data\retinograms\drive\test\vessel_masks\'];
fov_path = [aroot, 'data\retinograms\drive\test\foveal_masks\'];

img_dir = dir([img_path,'*.png']);
staal_dir = dir([staal_path,'*.png']);
msk_dir = dir([msk_path,'*.mat']);
fov_dir = dir([fov_path,'*.mat']);

raw_values = [];
staal_values = [];
classes = [];

figure(1); clf; hold on;
figure(2); clf; hold on;
for i = 1:5
    raw_img = double(imread([img_path, img_dir(i).name]));
    raw_img = 255-raw_img(:,:,2);
    raw_img = raw_img - min(raw_img(:));
    raw_img = raw_img / max(raw_img(:));
    
    staal_img = mean(imread([staal_path, staal_dir(i).name]),3);
    v_mask = u_load([msk_path, msk_dir(i).name]);
    f_mask = u_load([fov_path, fov_dir(i).name]);
    v_mask = v_mask & f_mask;
    
    raw_values = [raw_values; raw_img(f_mask)];
    staal_values = [staal_values; staal_img(f_mask)];
    classes = [classes; v_mask(f_mask)];

    [raw_roc, raw_Az] = ...
        calculate_roc_curve(raw_img(f_mask), v_mask(f_mask));
    figure(1);
        plot(raw_roc(:,1),raw_roc(:,2),'b-');

    [staal_roc, staal_Az] = ...
        calculate_roc_curve(staal_img(f_mask), v_mask(f_mask));
    figure(2);
        plot(staal_roc(:,1),staal_roc(:,2),'b-');
    
    disp([raw_Az staal_Az]);
end
    
lw = 2;

[raw_roc, raw_Az] = calculate_roc_curve(raw_values, classes);
figure(1);
    plot(raw_roc(:,1),raw_roc(:,2),'r-', 'linewidth', lw);
    axis([0,1,0,1],'square');
    
[staal_roc, staal_Az] = calculate_roc_curve(staal_values, classes);
figure(2);
    plot(staal_roc(:,1),staal_roc(:,2),'r-', 'linewidth', lw);
    axis([0,1,0,1],'square');
    
disp([raw_Az staal_Az]);


    
