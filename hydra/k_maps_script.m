function k_maps_script(idx)

%Script to compute results for k_maps expereiments on hydra

if ~ispc
     [z num_jobs] = unix('echo $NUM_JOBS'); num_jobs = str2num(num_jobs);
     [z mammo_data] = unix('echo $IMAGE_DIR'); mammo_data(end) = [];
end

%image_dir = [asymmetryroot 'data/' user_dir];
image_dir = [asymmetryroot 'data/mammograms/2004_screening_processed/' mammo_data '/'];
image_list = dir([image_dir '*.mat']);

%------------------------------------------------
% Generic steps to workout which images to process in each job
num_images = length(image_list);
images_per_job = ceil(num_images / num_jobs);
start_idx = (idx-1)*images_per_job + 1;
end_idx = min(num_images, idx*images_per_job);

%------------------------------------------------
% Any job specific computing here...
mam_names = get_mammo_info(image_list);
px_per_mm = 100 / 9;

lines_dir_all = [asymmetryroot 'data/line_maps/g2d_all/2004_screening_processed/' mammo_data '/'];
ori_dir_all = [asymmetryroot 'data/orientation_maps/g2d_all/2004_screening_processed/' mammo_data '/'];
lines_dir_1 = [asymmetryroot 'data/line_maps/g2d_1/2004_screening_processed/' mammo_data '/'];
ori_dir_1 = [asymmetryroot 'data/orientation_maps/g2d_1/2004_screening_processed/' mammo_data '/'];
lines_dir_2 = [asymmetryroot 'data/line_maps/g2d_2/2004_screening_processed/' mammo_data '/'];
ori_dir_2 = [asymmetryroot 'data/orientation_maps/g2d_2/2004_screening_processed/' mammo_data '/'];
lines_dir_3 = [asymmetryroot 'data/line_maps/g2d_3/2004_screening_processed/' mammo_data '/'];
ori_dir_3 = [asymmetryroot 'data/orientation_maps/g2d_3/2004_screening_processed/' mammo_data '/'];

mkdir(lines_dir_all);
mkdir(ori_dir_all);
mkdir(lines_dir_1);
mkdir(ori_dir_1);
mkdir(lines_dir_2);
mkdir(ori_dir_2);
mkdir(lines_dir_3);
mkdir(ori_dir_3);


display(['Processing images ' num2str(start_idx) ' to ' num2str(end_idx)]);
for ii = start_idx:end_idx
    display(['processing image ' num2str(ii)]);
    
    %Process the images
    
    mammogram = u_load([image_dir image_list(ii).name]);
    
    %All scales
%     [line_orientation, line_map] = karssemeijer_line_detection(...
%         mammogram,...
%         'line_scales', px_per_mm*[0.1 0.17 0.29],...
%         'grad_scale', px_per_mm,...
%         'grad_ori_thresh', pi/6,...
%         'grad_strength_thresh', 25,...
%         'line_strength_thresh', 0,...
%         'binary_map', 1,...
%         'degrees', 0);
%         
%     save([lines_dir_all mam_names{ii} '_lines.mat'], 'line_map');
%     save_uint8([ori_dir_all mam_names{ii} '_ori.mat'], line_orientation);
%     
%     %sigma 1
%     [line_orientation, line_map] = karssemeijer_line_detection(...
%         mammogram,...
%         'line_scales', px_per_mm*0.1,...
%         'grad_scale', px_per_mm,...
%         'grad_ori_thresh', pi/6,...
%         'grad_strength_thresh', 25,...
%         'line_strength_thresh', 0,...
%         'binary_map', 1,...
%         'degrees', 0);
%         
%     save([lines_dir_1 mam_names{ii} '_lines.mat'], 'line_map');
%     save_uint8([ori_dir_1 mam_names{ii} '_ori.mat'], line_orientation);
%     
%     %sigma 2
%     [line_orientation, line_map] = karssemeijer_line_detection(...
%         mammogram,...
%         'line_scales', px_per_mm*0.17,...
%         'grad_scale', px_per_mm,...
%         'grad_ori_thresh', pi/6,...
%         'grad_strength_thresh', 25,...
%         'line_strength_thresh', 0,...
%         'binary_map', 1,...
%         'degrees', 0);
%         
%     save([lines_dir_2 mam_names{ii} '_lines.mat'], 'line_map');
%     save_uint8([ori_dir_2 mam_names{ii} '_ori.mat'], line_orientation);

    %sigma 3
    [line_orientation, line_map] = karssemeijer_line_detection(...
        mammogram,...
        'line_scales', px_per_mm*0.29,...
        'grad_scale', px_per_mm,...
        'grad_ori_thresh', pi/6,...
        'grad_strength_thresh', 25,...
        'line_strength_thresh', 0,...
        'binary_map', 1,...
        'degrees', 0);
        
    save([lines_dir_3 mam_names{ii} '_lines.mat'], 'line_map');
    save_uint8([ori_dir_3 mam_names{ii} '_ori.mat'], line_orientation);


end