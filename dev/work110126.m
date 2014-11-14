ori_errors = [];
for ii = 1:20
    
    load(['Z:\data\synthetic_lines\philtres\image' zerostr(ii,3) '.mat']);
    ori_map = load_uint8(['Z:\data\synthetic_lines\philtres\results\191934\image' zerostr(ii,3) '_class.mat']);
    ori_map(label ~= 1) = nan;
    ori_errors = [ori_errors;...
        abs(mb_mod(ori_map(label_centre & label == 1) - label_orientation(label_centre & label == 1), 180))];
    
    
%     figure; 
%     subplot(1,2,1); imagesc(label_orientation); colormap([0 0 0; hsv(180)]); axis image; caxis([0 180])
%     subplot(1,2,2); imagesc(ori_map); colormap([0 0 0; hsv(180)]); axis image; caxis([0 180])
%     xlabel(['Mean absolute error = ' num2str(mean(ori_errors))]);
    clear label*
end
%%
ori_errors2 = [];
for ii = 1:100
    
    load(['C:\isbe\asymmetry_project\data\synthetic_lines\lines512\labels\label' zerostr(ii,3) '.mat']);
    ori_map = load_uint8(['C:\isbe\asymmetry_project\data\synthetic_lines\lines512\results\191934\probability_image' zerostr(ii,3) '.mat']);
    ori_map(label ~= 1) = nan;
    ori_errors2 = [ori_errors2;...
        abs(mb_mod(ori_map(label_centre & label == 1) - label_orientation(label_centre & label == 1), 180))];
    
    
%     figure; 
%     subplot(1,2,1); imagesc(label_orientation); colormap([0 0 0; hsv(180)]); axis image; caxis([0 180])
%     subplot(1,2,2); imagesc(ori_map); colormap([0 0 0; hsv(180)]); axis image; caxis([0 180])
%     xlabel(['Mean absolute error = ' num2str(mean(ori_errors))]);
    clear label*
end
%%
ori_errors_un = [];
ori_errors_no = [];
for ii = 1:20
    
    load(['Z:\data\synthetic_lines\philtres\image' zerostr(ii,3) '.mat']);
    ori_map = load_uint8(['Z:\data\synthetic_lines\philtres\results\233908\image' zerostr(ii,3) '_class.mat']);
    ori_map(label ~= 1) = nan;
    ori_errors_un = [ori_errors_un;...
        abs(mb_mod(ori_map(label_centre & label == 1) - label_orientation(label_centre & label == 1), 180))];
    
    load(['Z:\data\synthetic_lines\philtres\image' zerostr(ii,3) '.mat']);
    ori_map = load_uint8(['Z:\data\synthetic_lines\philtres\results\238797\image' zerostr(ii,3) '_class.mat']);
    ori_map(label ~= 1) = nan;
    ori_errors_no = [ori_errors_no;...
        abs(mb_mod(ori_map(label_centre & label == 1) - label_orientation(label_centre & label == 1), 180))];
    
%     figure; 
%     subplot(1,2,1); imagesc(label_orientation); colormap([0 0 0; hsv(180)]); axis image; caxis([0 180])
%     subplot(1,2,2); imagesc(ori_map); colormap([0 0 0; hsv(180)]); axis image; caxis([0 180])
%     xlabel(['Mean absolute error = ' num2str(mean(ori_errors3))]);
    clear label*
end
%%
%%
ori_errors5 = [];
for ii = 1:20
    
    load(['C:\isbe\asymmetry_project\data\synthetic_lines\curves512\labels\label' zerostr(ii,3) '.mat']);
    ori_map = load_uint8(['Z:\data\synthetic_lines\curves512\results\191934\image' zerostr(ii,3) '_class.mat']);
    ori_map(label ~= 1) = nan;
    ori_errors5 = [ori_errors5;...
        abs(mb_mod(ori_map(label_centre & label == 1) - label_orientation(label_centre & label == 1), 180))];
    
%     figure; 
%     subplot(1,2,1); imagesc(label_orientation); colormap([0 0 0; hsv(180)]); axis image; caxis([0 180])
%     subplot(1,2,2); imagesc(ori_map); colormap([0 0 0; hsv(180)]); axis image; caxis([0 180])
%     xlabel(['Mean absolute error = ' num2str(mean(ori_errors3))]);
    clear label*
end
%%
ori_errors4 = [];
for ii = 1:100
    
    load(['C:\isbe\asymmetry_project\data\synthetic_lines\real512\labels\label' zerostr(ii,3) '.mat']);
    ori_map = load_uint8(['C:\isbe\asymmetry_project\data\synthetic_lines\real512\results\233908\image' zerostr(ii,3) '_class.mat']);
    ori_map(label ~= 1) = nan;
    ori_errors4 = [ori_errors4;...
        abs(mb_mod(ori_map(label_centre & label == 1) - label_orientation(label_centre & label == 1), 180))];
    
    
%     figure; 
%     subplot(1,2,1); imagesc(label_orientation); colormap([0 0 0; hsv(180)]); axis image; caxis([0 180])
%     subplot(1,2,2); imagesc(ori_map); colormap([0 0 0; hsv(180)]); axis image; caxis([0 180])
%     xlabel(['Mean absolute error = ' num2str(mean(ori_errors))]);
    clear label*
end
%%
par_phil = u_load('Z:\data\synthetic_lines\philtres\parameters.mat');
angles = [];
widths = [];
contrasts = [];
squashes = [];
for ii = 1:100
    for jj = 1:length(par_phil(ii).curr_para)
            
        angles(end+1,1) = mod(par_phil(ii).curr_para(jj).orientation,180);
        widths(end+1,1) = par_phil(ii).curr_para(jj).halfwidth;
        contrasts(end+1,1) = par_phil(ii).curr_para(jj).contrast;
        squashes(end+1,1) = par_phil(ii).curr_para(jj).squash;
    end
end
%
par_mike = u_load('Z:\data\synthetic_lines\lines512\parameters.mat');
angles2 = [];
widths2 = [];
contrasts2 = [];
squashes2 = [];
for ii = 1:100
    for jj = 1:length(par_mike(ii).curr_para)
            
        angles2(end+1,1) = mod(par_mike(ii).curr_para(jj).orientation,180);
        widths2(end+1,1) = par_mike(ii).curr_para(jj).halfwidth;
        contrasts2(end+1,1) = par_mike(ii).curr_para(jj).contrast;
        squashes2(end+1,1) = par_mike(ii).curr_para(jj).squash;
    end
end