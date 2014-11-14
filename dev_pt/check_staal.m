% check Staal's Az

clear;

aroot = asymmetryroot;
img_path = [aroot, 'data\retinograms\drive\test\output\Staal_etal\'];
msk_path = [aroot, 'data\retinograms\drive\test\vessel_masks\'];
fov_path = [aroot, 'data\retinograms\drive\test\foveal_masks\'];

img_dir = dir([img_path,'*.png']);
msk_dir = dir([msk_path,'*.mat']);
fov_dir = dir([fov_path,'*.mat']);

values = [];
classes = [];

for i = 1:20
    img = mean(imread([img_path, img_dir(i).name]),3);
    v_mask = u_load([msk_path, msk_dir(i).name]);
    f_mask = u_load([fov_path, fov_dir(i).name]);
    v_mask = v_mask & f_mask;
    
    values = [values; img(f_mask)];
    classes = [classes; v_mask(f_mask)];
end
    
[roc, Az] = calculate_roc_curve(values, classes);

figure(10); clf; hold on;
    plot(roc(:,1),roc(:,2),'b-');
%     plot(roc2(:,1),roc2(:,2),'r-');
    axis([0,1,0,1],'square');
Az


    
