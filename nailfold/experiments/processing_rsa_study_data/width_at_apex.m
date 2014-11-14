function [] = width_at_apex(start_i, end_i)
%Looking at remeasuring apex width

if ~exist('start_i', 'var') || isempty(start_i)
    start_i = 1;
end
if ~exist('end_i', 'var') || isempty(end_i)
    end_i = start_i;
end

args.data_dir = [nailfoldroot 'data/rsa_study/test_half/'];
args.image_dir = 'images';
args.prob_dir = 'rf_classification/296655';
args.ori_dir = 'rf_regression/296621';
args.width_dir = 'rf_regression/297037';
args.candidates_dir = 'apex_maps\set12g_half_296655\mixed_maxima_new\selected_candidates';
args.apex_gt_dir = 'apex_gt';
args.prob_sigma = 2;
args.ori_sigma = 0;
args.width_sigma = 2;
args.base_width = 20;

image_dir = [args.data_dir args.image_dir '/'];
prob_dir = [args.data_dir 'predictions/detection/' args.prob_dir '/'];
ori_dir = [args.data_dir 'predictions/orientation/' args.ori_dir '/'];
width_dir = [args.data_dir 'predictions/width/' args.width_dir '/'];
candidates_dir = [args.data_dir '/' args.candidates_dir '/'];
apex_gt_dir = [args.data_dir '/' args.apex_gt_dir '/'];

pred_list = dir([prob_dir '*.mat']);
num_images = length(pred_list);

% %1. Workout number of images in job
% job_size = ceil(num_images / args.num_jobs);
% 
% %2. Workout start and end indices for job
% start_i	= (args.task_id-1)*job_size + 1;
% end_i	= min(args.task_id*job_size, num_images);
display(['Extracting vessel centres from images ' num2str(start_i) ' to ' num2str(end_i)]);

x = repmat(-31.5:31.5, 128, 1);
y = repmat((-31.5:95.5)', 1, 64);
xy = [x(:) y(:)];
patch_sz = size(x);

%Loop though each image
for i_im = start_i:end_i
    im_name = pred_list(i_im).name(1:end-9);
    
    display(['Processing image ' num2str(i_im) ', ' datestr(now)]);  
    
    %Load in images and vessel markup
    vessel_im = u_load([image_dir im_name '.mat']);
    vessel_prob = u_load([prob_dir pred_list(i_im).name]);
    vessel_ori = u_load([ori_dir pred_list(i_im).name]);
    vessel_width = u_load([width_dir pred_list(i_im).name]);
    load([candidates_dir im_name '_candidates'],...
        'candidate_xy', 'kept'); 
    load([apex_gt_dir im_name '_gt.mat'], 'apex_xy', 'apex_widths', 'is_distal');
    
    candidate_xy = candidate_xy(kept,:);
    num_apexes = size(candidate_xy, 1);        
    
    [~, ~, candidate_detections] =...
        evaluate_apex_candidates(apex_xy, candidate_xy, apex_widths/2, vessel_prob, ...
        [], [], 0);
        
    %Smooth the vessel probs
    if args.prob_sigma
        g_prob = gaussian_filters_1d(args.prob_sigma);
        g_prob = g_prob / sum(g_prob);    
        vessel_prob = conv2(g_prob', g_prob, vessel_prob, 'same');
    end
    if args.ori_sigma
        g_ori = gaussian_filters_1d(args.ori_sigma);
        g_ori = g_ori / sum(g_ori);
        vessel_ori = conv2(g_ori', g_ori, vessel_ori, 'same');
    end 
    if args.width_sigma
        g_width = gaussian_filters_1d(args.width_sigma);
        g_width = g_width / sum(g_width);
        vessel_width = conv2(g_width', g_width, vessel_width, 'same');
    end   

    for i_pt = 1:min(num_apexes,20)
        
        if ~candidate_detections(i_pt) || ~is_distal(candidate_detections(i_pt))
            continue;
        end

        %Get predicted scale and orientation at this point
        ax = candidate_xy(i_pt,1);
        ay = candidate_xy(i_pt,2);
        a_ori = vessel_ori(round(ay), round(ax));
        a_theta = angle(a_ori) / 2;
        a_width = vessel_width(round(ay), round(ax));
        apex_measures.base_orientation(i_pt) = a_ori;
    
        %Get scale relative to base width a make rotation matrix
        rot = [cos(a_theta) -sin(a_theta); sin(a_theta) cos(a_theta)];
        scale = a_width / args.base_width;

        %Transform points given scale and angle and translate to
        %candidate position
        xya = xy * rot * scale;
        xa = reshape(xya(:,1) + ax, patch_sz);
        ya = reshape(xya(:,2) + ay, patch_sz);
        
        gt_width = apex_widths(candidate_detections(i_pt)) / (2*scale);
        a_width = a_width / scale;
        
        %Sample vessel prob patch
        vessel_im_patch = interp2(vessel_im, xa, ya, '*linear', 0);
        vessel_prob_patch = interp2(vessel_prob, xa, ya, '*linear', 0);
        vessel_ori_patch = interp2(vessel_ori, xa, ya, '*linear', 0);
        vessel_width_patch = interp2(vessel_width, xa, ya, '*linear', 0);
        vessel_ori_patch = vessel_ori_patch * conj(a_ori);
        
        figure;
        subplot(1,3,1);
        imgray(vessel_im_patch);
        plot([patch_sz(2) patch_sz(2)]/2 + 1, (patch_sz(2)+[-1 1]*gt_width)/2, 'r-x');
        plot([patch_sz(2) patch_sz(2)]/2 - 1, (patch_sz(2)+[-1 1]*a_width)/2, 'g-x');
        
        subplot(1,3,2);
        imgray(vessel_prob_patch);
        plot([patch_sz(2) patch_sz(2)]/2 + 1, (patch_sz(2)+[-1 1]*gt_width)/2, 'r-x');
        plot([patch_sz(2) patch_sz(2)]/2 - 1, (patch_sz(2)+[-1 1]*a_width)/2, 'g-x');
        
        subplot(1,3,3);
        imgray(complex2rgb(vessel_ori_patch));
    end
end