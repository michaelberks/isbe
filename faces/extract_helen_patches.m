function [] = extract_helen_patches(im_name, anno_dir, image_dir, output_dir, list_fid, do_plot)
%EXTRACT_HELEN_PATCHES *Insert a one line summary here*
%   [] = extract_helen_patches()
%
% Inputs:
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 25-May-2016
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
if ~exist('output_dir', 'var')
    output_dir = [];
end
if ~exist('list_fid', 'var')
    list_fid = [];
end
if ~exist('do_plot', 'var')
    do_plot = 0;
end

%Define indices to face objects - bit annoying to repeat for each image but
%only takes a second TO DO set up uparkargs list
x = repmat(-63.5:63.5, 128, 1);
y = x';
xy = [x(:) y(:)];

base_iod = 100;

eye_l_idx = 135:154;
eye_r_idx = 115:134;
nose_idx = 42:58;
mouth_idx = 59:87;

%Load face and annotation
face_pts = load([anno_dir im_name '.txt']);
face_im = double(imread([image_dir im_name '.jpg']));

eye_l_centre = mean(face_pts(eye_l_idx,:));
eye_r_centre = mean(face_pts(eye_r_idx,:));
nose_centre = mean(face_pts(nose_idx,:));
mouth_centre = mean(face_pts(mouth_idx,:));

%Compute transform matrix based on inter-occular distance
inter_oc_vec = diff(face_pts([135 115],:));

t_mat = [inter_oc_vec; -inter_oc_vec(2), inter_oc_vec(1)] / base_iod;
xyi = xy*t_mat;

%Extract patches from left and right eyes, nose and mouth
eye_l_patch = uint8(cat(3,...
    reshape(interp2(face_im(:,:,1), xyi(:,1) + eye_l_centre(1), xyi(:,2) + eye_l_centre(2), '*linear'), 128, 128),...
    reshape(interp2(face_im(:,:,2), xyi(:,1) + eye_l_centre(1), xyi(:,2) + eye_l_centre(2), '*linear'), 128, 128),...
    reshape(interp2(face_im(:,:,3), xyi(:,1) + eye_l_centre(1), xyi(:,2) + eye_l_centre(2), '*linear'), 128, 128)));

eye_r_patch = uint8(cat(3,...
    reshape(interp2(face_im(:,:,1), xyi(:,1) + eye_r_centre(1), xyi(:,2) + eye_r_centre(2), '*linear'), 128, 128),...
    reshape(interp2(face_im(:,:,2), xyi(:,1) + eye_r_centre(1), xyi(:,2) + eye_r_centre(2), '*linear'), 128, 128),...
    reshape(interp2(face_im(:,:,3), xyi(:,1) + eye_r_centre(1), xyi(:,2) + eye_r_centre(2), '*linear'), 128, 128)));

nose_patch = uint8(cat(3,...
    reshape(interp2(face_im(:,:,1), xyi(:,1) + nose_centre(1), xyi(:,2) + nose_centre(2), '*linear'), 128, 128),...
    reshape(interp2(face_im(:,:,2), xyi(:,1) + nose_centre(1), xyi(:,2) + nose_centre(2), '*linear'), 128, 128),...
    reshape(interp2(face_im(:,:,3), xyi(:,1) + nose_centre(1), xyi(:,2) + nose_centre(2), '*linear'), 128, 128)));

mouth_patch = uint8(cat(3,...
    reshape(interp2(face_im(:,:,1), xyi(:,1) + mouth_centre(1), xyi(:,2) + mouth_centre(2), '*linear'), 128, 128),...
    reshape(interp2(face_im(:,:,2), xyi(:,1) + mouth_centre(1), xyi(:,2) + mouth_centre(2), '*linear'), 128, 128),...
    reshape(interp2(face_im(:,:,3), xyi(:,1) + mouth_centre(1), xyi(:,2) + mouth_centre(2), '*linear'), 128, 128)));
  
%Save
if ~isempty(output_dir)
    imwrite(eye_l_patch, [output_dir im_name '_el.png']);
    imwrite(eye_r_patch, [output_dir im_name '_er.png']);
    imwrite(nose_patch, [output_dir im_name '_no.png']);
    imwrite(mouth_patch, [output_dir im_name '_mo.png']);
end

%List
if ~isempty(list_fid)
    fprintf(list_fid, '%s %d \n', [im_name '_el.png'], 0);
    fprintf(list_fid, '%s %d \n', [im_name '_er.png'], 1);
    fprintf(list_fid, '%s %d \n', [im_name '_no.png'], 2);
    fprintf(list_fid, '%s %d \n', [im_name '_mo.png'], 3);
end

%Plot
if do_plot
    figure;
    subplot(2,2,1); imgray(eye_l_patch);
    subplot(2,2,2); imgray(eye_r_patch);
    subplot(2,2,3); imgray(nose_patch);
    subplot(2,2,4); imgray(mouth_patch);
end
