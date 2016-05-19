function [] = make_apex_candidate_ims(varargin)
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
args = u_packargs(varargin, 0, ... % the user's input
    {'image_names'},           ...
    'data_dir',             [nailfoldroot 'data/rsa_study/master_set/'],...
    'feature_im_dir',       'images',...
    'feature_im_ext',       '.mat',...
    'ori_dir',              'rf_regression/296621/',...
    'width_dir',            'rf_regression/297037/',...
    'feature_sigma',        0,...
    'ori_sigma',            0,...
    'width_sigma',          2,...
    'candidates_dir',       'apex_maps/local_maxima',...
    'output_im_dir',        'apex_images',...
    'output_format',        'graysc',...
    'patch_lims',           [-63.5 63.5 -31.5 95.5],...
    'base_width',           20, ...
    'overwrite',            0,...
    'debug',                0);

%Form full directory paths and create folder for HoGs
feature_im_dir = [args.data_dir args.feature_im_dir '/'];
ori_dir = [args.data_dir 'predictions/orientation/' args.ori_dir '/'];
width_dir = [args.data_dir 'predictions/width/' args.width_dir '/'];
candidates_dir = [args.data_dir '/' args.candidates_dir '/'];
output_im_dir = [args.data_dir '/' args.output_im_dir '/'];
create_folder(output_im_dir);

%Set up x,y coordinates for template patch
x = args.patch_lims(1):args.patch_lims(2);
y = args.patch_lims(3):args.patch_lims(4);
size_x = length(x);
size_y = length(y);
xx = repmat(x, size_y, 1);
yy = repmat(y', 1, size_x);
xy = [xx(:) yy(:)];

if args.feature_sigma
    g_feat = gaussian_filters_1d(args.feature_sigma);
    g_feat = g_feat / sum(g_feat);    
end
if args.ori_sigma
    g_ori = gaussian_filters_1d(args.ori_sigma);
    g_ori = g_ori / sum(g_ori);
end 
if args.width_sigma
    g_width = gaussian_filters_1d(args.width_sigma);
    g_width = g_width / sum(g_width);
end 

num_images = length(args.image_names);
%Loop though each image
for i_im = 1:num_images
    
    im_name = args.image_names{i_im} ;  
    display(['Processing image ' num2str(i_im) ', ' datestr(now)]);  
    
    %Load in images and vessel markup
    vessel_feature_im = double(u_load([feature_im_dir im_name args.feature_im_ext]));
    vessel_ori = u_load([ori_dir im_name '_pred.mat']);
    vessel_width = u_load([width_dir im_name '_pred.mat']);
    load([candidates_dir im_name '_candidates'], 'candidate_xy');
        
    if isempty(candidate_xy) %#ok
        candidates_hogs = []; %#ok
        candidate_oris = []; %#ok
        candidate_widths = []; %#ok
    else
    
        if args.feature_sigma
            vessel_feature_im = conv2(g_feat', g_feat, vessel_feature_im, 'same');
        end
        if args.ori_sigma
            vessel_ori = conv2(g_ori', g_ori, vessel_ori, 'same');
        end 
        if args.width_sigma
            vessel_width = conv2(g_width', g_width, vessel_width, 'same');
        end
        candidate_oris = interp2(vessel_ori, candidate_xy(:,1), candidate_xy(:,2));
        candidate_widths = interp2(vessel_width, candidate_xy(:,1), candidate_xy(:,2));
        clear vessel_ori vessel_width;

        candidate_oris = angle(candidate_oris) / 2; 
        num_candidates = size(candidate_xy,1);            

        for i_can = 1:num_candidates

            %EXtract patch and compute HoG
            ori_c = candidate_oris(i_can);
            width_c = candidate_widths(i_can);

            %Get scale relative to base width a make rotation matrix
            rot = [cos(ori_c) -sin(ori_c); sin(ori_c) cos(ori_c)];
            scale = width_c / args.base_width;

            %Transform points given scale and angle and translate to
            %candidate position
            xya = xy * rot * scale;
            xa = reshape(xya(:,1) + candidate_xy(i_can,1), size_y, size_x);
            ya = reshape(xya(:,2) + candidate_xy(i_can,2), size_y, size_x);

            %Sample vessel prob patch
            vessel_feature_patch = interp2(vessel_feature_im, xa, ya, '*linear', 0);
            
            patch_name = [output_im_dir im_name '_' zerostr(i_can,3) '.png'];
            
            if ~args.overwrite && exist(patch_name, 'file')
                display(['Skipping image ' num2str(i_im) ',' patch_name ' already exists. Re-run with overwrite=1 if necessary']); 
                continue;
            end
    
            switch args.output_format
                case 'gray'
                    imwrite(vessel_feature_patch, patch_name);
                case 'graysc'
                    write_im_from_colormap(vessel_feature_patch, patch_name, gray(256));
                case 'complex'
                    imwrite(complex2rgb(vessel_feature_patch), patch_name);
                otherwise
                    display(['Invalid ouput format: ' args.output_format ', should be gray, graysc or complex']);
            end
                    
        end
    end
end          
    
    
