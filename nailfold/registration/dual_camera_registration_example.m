%Example script for preparing and the registering frames captured from the
%dual camera system

%Set folder containing all the frames
frame_set = {...
    'C:\isbe\nailfold\data\dual_wavelength\A Tonia 3 secs occ 1\2013_09_12\';...
    'C:\isbe\nailfold\data\dual_wavelength\A Tonia 3 secs rel 1\2013_09_12\';...
    'C:\isbe\nailfold\data\dual_wavelength\A Tonia 3 secs baseline 1\2013_09_12\';...
    'C:\isbe\nailfold\data\dual_wavelength\A Tonia 3 secs post 1\2013_09_12\';...
    'C:\isbe\nailfold\data\dual_wavelength\2013_07_12\original\' };
frames_dir = frame_set{4};

%Prepare the frames by removing the dirt from each camera and computing how
%all the consecutive frames within each camera register to each other -
%warning, this may take a while for long sequences.
prepare_dual_camera_frames(frames_dir,...
    'rotation1', 1,...
    'rotation2', 1,...
    'flip1', 0,...
    'flip2', 1);

%Load in the transforms from the above registrations
c1_transforms = u_load([frames_dir 'transforms\camera1_final_reg_transforms.mat']);
c2_transforms = u_load([frames_dir 'transforms\camera2_final_reg_transforms.mat']);

%Register pairs of compund frames from each camera - define the range used
%in each compound frame either by frame numbers (which will be the same for
%both cameras) or by time interval in seconds over the sequence (in which
%case different frames may be used for each camera because the times may
%not match up)
register_dual_camera_frames([frames_dir 'final_corrections\'],...
    c1_transforms, c2_transforms, [1 30], 'frames',...
    'save_images', 1, ...,
    'save_dir', [frames_dir 'difference_images_f\']);

register_dual_camera_frames([frames_dir 'final_corrections\'],...
    c1_transforms, c2_transforms, [0 1], 'time',...
    'save_images', 1, ...
    'save_dir', [frames_dir 'difference_images_t\']);
%%
load([frames_dir 'difference_images_f\difference_image001.mat']);

rf_det = u_load('C:\isbe\nailfold\models\vessel\detection\rf_classification\296655\predictor.mat');
rf_det.tree_root = 'C:\isbe\nailfold\models\vessel\detection\rf_classification/';
det_args = u_load('C:\isbe\nailfold\models\vessel\detection\rf_classification\296655\job_args.mat');

valid_mask = ~isnan(registered_mosaics(:,:,1)) & ~isnan(registered_mosaics(:,:,2));

vessel_im1 = registered_mosaics(:,:,1);
vessel_im1(isnan(vessel_im1)) = mean(vessel_im1(valid_mask));
vessel_im2 = registered_mosaics(:,:,2);
vessel_im2(isnan(vessel_im2)) = mean(vessel_im2(valid_mask));

[vessel_prob1] = predict_image(...
    'image_in',vessel_im1,...
    'mask', imerode(valid_mask, strel('disk', 10)), ...
    'decomposition_args', det_args.decomposition_args,...
    'predictor', rf_det, ...
    'prediction_type', 'rf_classification',...
    'output_type', 'detection');

[vessel_prob2] = predict_image(...
    'image_in',vessel_im2,...
    'mask', imerode(valid_mask, strel('disk', 10)), ...
    'decomposition_args', det_args.decomposition_args,...
    'predictor', rf_det, ...
    'prediction_type', 'rf_classification',...
    'output_type', 'detection');

vessel_mask = vessel_prob1 > 0.8 & vessel_prob2 > 0.8;

figure; 
subplot(1,2,1); imgray(vessel_prob1);
subplot(1,2,2); imgray(vessel_prob2);


save([frames_dir 'difference_images_f\vessel_masks001.mat'], 'vessel_prob1', 'vessel_prob2', 'vessel_mask'); 
%%
frame_set = {...
    'C:\isbe\nailfold\data\dual_wavelength\A Tonia 3 secs occ 1\2013_09_12\';...
    'C:\isbe\nailfold\data\dual_wavelength\A Tonia 3 secs rel 1\2013_09_12\';...
    'C:\isbe\nailfold\data\dual_wavelength\A Tonia 3 secs baseline 1\2013_09_12\';...
    'C:\isbe\nailfold\data\dual_wavelength\A Tonia 3 secs post 1\2013_09_12\';...
    'C:\isbe\nailfold\data\dual_wavelength\2013_07_12\original\' };

for i_fr = 1:4
    frames_dir = frame_set{i_fr};

    load([frames_dir 'difference_images_f\vessel_masks001.mat'], 'vessel_prob1', 'vessel_prob2', 'vessel_mask'); 
    load([frames_dir 'difference_images_f\difference_image001.mat']);
    
    valid_mask = ~isnan(registered_mosaics(:,:,1)) & ~isnan(registered_mosaics(:,:,2));
    %
    bias = 0;%mean(registered_difference(valid_mask & ~vessel_mask));
    vessel_difference = registered_difference;
    vessel_difference(~vessel_mask) = bias;
    
    [pos_y pos_x] = find(vessel_difference > bias + 1);
    [neg_y neg_x] = find(vessel_difference < bias - 1);
    
    vessel_im = registered_difference;
    vessel_im(~vessel_mask) = -inf;    
    
    cmap = [1 1 1; jet(255)];
    figure; imgray(vessel_mask);
    figure; imgray(registered_difference);  caxis([-10 10]); colorbar;%caxis([min(registered_difference(valid_mask)) max(registered_difference(valid_mask))]);
    figure; imgray(vessel_im); colormap(cmap); caxis([-10 10]); colorbar;

    figure; imgray(mean(registered_mosaics,3));
    plot(pos_x, pos_y, 'g.');
    plot(neg_x, neg_y, 'r.');
end
