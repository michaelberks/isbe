function [] = compute_vessel_centre_hogs(varargin)
%EXTRACT_VESSEL_CENTRES *Insert a one line summary here*
%   [] = extract_vessel_centres()
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
% Created: 18-Jun-2013
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
args = u_packargs(varargin,... % the user's input
    '0', ... % non-strict mode
    'task_id',              unixenv('SGE_TASK_ID',1), ...
    'num_jobs',             unixenv('NUM_JOBS',100), ...
    'data_dir',             [nailfoldroot 'data/rsa_study/set12g/'],...
    'feature_im_dir',       'predictions/detection/rf_classification/257273/',...
    'centre_dir',           'vessel_centres/',...
    'hog_dir',              'vessel_hogs/',...
    'max_size',             unixenv('MAX_SIZE', 1000),...
    'smoothing_sigma',      2,...
    'num_cells',            8,...
    'cell_sz',              8,... %Size of HoG cells in blocks
    'block_sz',             [2 2],...%Size of blocks in cells
    'num_ori_bins',         9,... %Number of bins in orientation histograms
    'norm_method',          'l1-sqrt',... %Method for local normalisation
    'block_spacing',        8,... %Separation between blocks in pixels - controls level of overlap
    'gradient_operator',    [-1 0 1],...
    'spatial_sigma',        0, ...
    'angle_wrap',           1,...
    'base_width',           20, ...
    'dist_thresh',          24^2,...
    'make_parts_folders',   0);

%Form full directory paths and create folder for HoGs
feature_im_dir = [args.data_dir args.feature_im_dir];
centre_dir = [args.data_dir args.centre_dir];
hog_dir = [args.data_dir args.hog_dir ];
create_folder(hog_dir);

%Get HoG args from main args
hog_args.cell_sz = [args.cell_sz args.cell_sz];
hog_args.block_sz = args.block_sz;
hog_args.num_ori_bins = args.num_ori_bins;
hog_args.norm_method = args.norm_method;
hog_args.block_spacing = args.block_spacing;
hog_args.gradient_operator = args.gradient_operator;
hog_args.spatial_sigma = args.spatial_sigma;
hog_args.angle_wrap = args.angle_wrap;

%1. Workout number of images in job
image_list = dir([feature_im_dir '*.mat']);
num_images = length(image_list);
job_size = ceil(num_images / args.num_jobs);

%2. Workout start and end indices for job
start_i	= (args.task_id-1)*job_size + 1;
end_i	= min(args.task_id*job_size, num_images);
display(['Computing HoGs for images ' num2str(start_i) ' to ' num2str(end_i)]);

%Create smoothing kernel for feature image
if args.smoothing_sigma
    g = gaussian_filters_1d(args.smoothing_sigma);
    g = g / sum(g);
end

%Get patch size and form template x,y coordinates for the patch
patch_sz = args.num_cells*args.cell_sz;
patch_sz = patch_sz + 2; %Account for padding
patch_sz2 = (patch_sz - 1)/2;
hog_sz = args.num_cells*args.num_cells*args.num_ori_bins;

%Set up x,y coordinates for template patch
x = repmat(-patch_sz2:patch_sz2, patch_sz, 1);
y = x';
xy = [x(:) y(:)];

%Loop though each image
for i_im = start_i:end_i
    im_name = image_list(i_im).name(1:end-9);
    
    display(['Processing image ' num2str(i_im) ', ' datestr(now)]);  
    
    %Load in images and vessel markup
    try
        vessel_feature_im = u_load([feature_im_dir image_list(i_im).name]);

        load([centre_dir im_name '_vc.mat'], 'vessel_centre', 'nrows', 'ncols');
        
    catch last_err
        display(last_err.message);
        continue;
    end
    
    %Smooth the vessel probs
    if args.smoothing_sigma
        vessel_feature_im = conv2(g', g, vessel_feature_im, 'same');
    end

    num_pts = length(vessel_centre.x);  

    if args.make_parts_folders
        %Make a folder for all the parts
        create_folder([hog_dir im_name '_hog']);
    end
    curr_pt = 1;
    part_num = 1;
    for i_pt = 1:num_pts
        
        %If the first point in a new part, make a container for the HoGs
        if curr_pt == 1
            part_sz = min(args.max_size, num_pts-i_pt+1);
            vessel_hog = zeros(part_sz, hog_sz);
        end

        %Get predicted scale and orientation at this point
        vxc = vessel_centre.x(i_pt);
        vyc = vessel_centre.y(i_pt);
        ori_c = angle(vessel_centre.ori(i_pt))/2;
        width_c = vessel_centre.width(i_pt);

        %Get scale relative to base width a make rotation matrix
        rot = [cos(ori_c) -sin(ori_c); sin(ori_c) cos(ori_c)];
        scale = width_c / args.base_width;

        %Transform points given scale and angle and translate to
        %candidate position
        xya = xy * rot * scale;
        xa = reshape(xya(:,1) + vxc, patch_sz, patch_sz);
        ya = reshape(xya(:,2) + vyc, patch_sz, patch_sz);

        %Sample vessel prob patch
        vessel_feature_patch = interp2(vessel_feature_im, xa, ya, '*linear', 0);
        [hog] = compute_HoG(vessel_feature_patch, hog_args);       
        vessel_hog(curr_pt,:) = hog(:)';
        
        %Increment the curr_pt and check if we've filled the part. If so, set
        %the curr_pt back to 1, save the HoGs, and increment the part number
        if curr_pt == part_sz
            if args.make_parts_folders                
                hog_name = [hog_dir im_name '_hog/'];
            else
                hog_name = [hog_dir im_name '_hog_'];                
            end
            save([hog_name 'part_' zerostr(part_num,4) '.mat'], 'vessel_hog');
            part_num = part_num + 1;
        else
            curr_pt = curr_pt + 1;
        end
    end       
end          
    
    
