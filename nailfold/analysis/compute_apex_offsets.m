function [] = compute_apex_offsets(varargin)

args = u_packargs(varargin,... % the user's input
    '0', ... % non-strict mode
    'task_id',              unixenv('SGE_TASK_ID',1), ...
    'num_jobs',             unixenv('NUM_JOBS',100), ...
    'data_dir',             [nailfoldroot 'data/rsa_study/set12g/'],...
    'centre_dir',           'vessel_centres/',...
    'contour_dir',          'vessel_contours/',...
    'mask_dir',             'vessel_masks/',...
    'apex_class_dir',       'apex_class_data/',...
    'base_width',           20,...
    'dist_thresh',          24);
    
%
centre_dir = [args.data_dir args.centre_dir];
centres_list = dir([centre_dir '*.mat']);
contour_dir = [args.data_dir args.contour_dir];
contours_list = dir([contour_dir '*.mat']);
mask_dir = [args.data_dir args.mask_dir];

apex_class_dir = [args.data_dir args.apex_class_dir];
create_folder(apex_class_dir);

num_ims = length(centres_list);
if length(contours_list) ~= num_ims
    error('Number of centres and number of contours do not match!');
end

dist_thresh = args.dist_thresh^2;

%1. Workout number of images in job
job_size = ceil(num_ims / args.num_jobs);

%2. Workout start and end indices for job
start_i	= (args.task_id-1)*job_size + 1;
end_i	= min(args.task_id*job_size, num_ims);
display(['Extracting vessel centres from images ' num2str(start_i) ' to ' num2str(end_i)]);

%Loop though each image
for i_im = start_i:end_i

    display(['Processing image ' num2str(i_im) ' of ' num2str(num_ims)]);

    %Load in the data we need
    s = load([centre_dir centres_list(i_im).name], 'vessel_centre');
    vessel_centre = s.vessel_centre; clear s;
    s = load([contour_dir contours_list(i_im).name], 'apex_idx', 'vessel_centre');
    vessel_contour = s.vessel_centre;
    apex_idx = s.apex_idx; clear s
    vessel_mask = u_load([mask_dir centres_list(i_im).name(1:end-7) '_v_mask.mat']);
    
    num_pts = length(vessel_centre.y);       

    %Correct apex coordinates frame
    apex_xy = [vessel_contour(apex_idx,1) vessel_contour(apex_idx,2)];

    apex_offsets = zeros(num_pts,2);
    apex_class = false(num_pts,1);

    for i_pt = 1:num_pts

        %Get predicted scale and orientation at this point
        vxc = vessel_centre.x(i_pt);
        vyc = vessel_centre.y(i_pt);
        ori_c = angle(vessel_centre.ori(i_pt))/2;
        width_c = vessel_centre.width(i_pt);

        %Get scale relative to base width a make rotation matrix
        rot = [cos(ori_c) -sin(ori_c); sin(ori_c) cos(ori_c)];
        scale = width_c / args.base_width;

        apex_dist = inf;
        for i_ap = 1:length(apex_idx)
            %Transform vessel centre
            apex_xy_t = (apex_xy(i_ap,:) - [vxc vyc])*rot'/scale;
            apex_dist_i = (sum(apex_xy_t.^2));

            if apex_dist_i < apex_dist
                apex_offsets(i_pt,:) = apex_xy_t;
                apex_dist = apex_dist_i;
            end
        end

        %Display patch if apex_xy_t lies within some dist
        if  vessel_mask(vyc, vxc) && (apex_dist < dist_thresh)  
            apex_class(i_pt) = 1;
        end
    end
    save([apex_class_dir centres_list(i_im).name(1:end-4) '_labels.mat'], 'apex_class', 'apex_offsets');

end