function [] = extract_helen_eye_patches(im_name, varargin)
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
args = u_packargs(varargin,... % the user's input
    '0', ... % non-strict mode
    'study_dir',             [toyotaroot 'data/helen'],...
    'annotation_dir',       'annotations',...
    'image_dir',            'images',...
    'output_dir',           [],...
    'list_fid',             [],...
    'n_samples_per_eye',    100,...
    'rand_sample_fun',      'uniform',...
    'rand_fun_params',      [100 100],...
    'patch_x',              -63.5:63.5,...
    'patch_y',              -63.5:63.5,...
    'do_left',              true,...
    'do_right',             true,...
    'mirror_right',         true,...
    'eye_l_idx',            135:154,...
	'eye_r_idx',            115:134,...
	'nose_idx',             42:58,...
	'mouth_idx',            59:87,...
    'iod_idx',              [135 115],...
    'base_iod',             100,...
	'plot',                 0);

%Set up folders
annotation_dir = [args.study_dir '/' args.annotation_dir '/'];
image_dir = [args.study_dir '/' args.image_dir '/'];

%Set up sampling points
[x y] = meshgrid(args.patch_x(:), args.patch_y(:));
xy = [x(:) y(:)];
patch_sz = size(x);

%Set up sampling function
switch(args.rand_sample_fun)
%                Criterion function   Is it an impurity measure?
%                ------------------   --------------------------
    case 'uniform',     samplefun = @uniform;
    case 'gaussian',   	samplefun = @gaussian;
    otherwise,     error([args.rand_sample_fun 'is not a recognised ''rand_sample_fun'' parameter.'])
end

%Load face and annotation
face_pts = load([annotation_dir im_name '.txt']);
face_im = double(imread([image_dir im_name '.jpg']));

%Compute transform matrix based on inter-occular distance
inter_oc_vec = diff(face_pts([135 115],:));

t_mat = [inter_oc_vec; -inter_oc_vec(2), inter_oc_vec(1)] / args.base_iod;
xyt = xy*t_mat;

xt = reshape(xyt(:,1), patch_sz);
yt = reshape(xyt(:,2), patch_sz);

patch_number = 1;
if args.do_left
    eye_l_centre = mean(face_pts(args.eye_l_idx,:));
    for i_pa = 1:args.n_samples_per_eye
        
        %sample offsets
        [ox, oy] = feval(samplefun, args.rand_fun_params);
        
        xi = xt + eye_l_centre(1) + ox;
        yi = yt + eye_l_centre(2) + oy;

        %Extract patches left eye
        eye_patch = uint8(cat(3,...
            interp2(face_im(:,:,1), xi, yi, '*linear'),...
            interp2(face_im(:,:,2), xi, yi, '*linear'),...
            interp2(face_im(:,:,3), xi, yi, '*linear')));
        
        patch_name = [im_name '_' zerostr(patch_number,4) '.png'];
        
        %Save
        if ~isempty(args.output_dir)
            imwrite(eye_patch, [args.output_dir patch_name]);
        end

        %List
        if ~isempty(args.list_fid)
            fprintf(args.list_fid, '%s %.2f %.2f \n', patch_name, ox, oy);
        end

        %Plot
        if args.plot >= i_pa
            figure; imgray(eye_patch);
            title('Left eye');
            xlabel(['x offset = ' num2str(ox, 3)]);
            ylabel(['y offset = ' num2str(oy, 3)]);
        end
        patch_number = patch_number + 1;
    end
end

if args.do_right
    eye_r_centre = mean(face_pts(args.eye_r_idx,:));
    
    %Flip coordinates horizontally to sample mirror image
    if args.mirror_right
        xt = fliplr(xt);
        yt = fliplr(yt);
    end
    
    for i_pa = 1:args.n_samples_per_eye
        
        %sample offsets
        [ox, oy] = feval(samplefun, args.rand_fun_params);
        
        xi = xt + eye_r_centre(1) + ox;
        yi = yt + eye_r_centre(2) + oy;

        %Extract patches left eye
        eye_patch = uint8(cat(3,...
            reshape(interp2(face_im(:,:,1), xi, yi, '*linear'), patch_sz),...
            reshape(interp2(face_im(:,:,2), xi, yi, '*linear'), patch_sz),...
            reshape(interp2(face_im(:,:,3), xi, yi, '*linear'), patch_sz)));
        
        if args.mirror_right
            %Record mirrored x-offset
            ox = -ox;
        end
        
        patch_name = [im_name '_' zerostr(patch_number,4) '.png'];
        
        %Save
        if ~isempty(args.output_dir)
            imwrite(eye_patch, [args.output_dir patch_name]);
        end

        %List
        if ~isempty(args.list_fid)
            fprintf(args.list_fid, '%s %.2f %.2f \n', patch_name, ox, oy);
        end

        %Plot
        if args.plot >= i_pa
            figure; imgray(eye_patch);
            title('Right eye');
            xlabel(['x offset = ' num2str(ox, 3)]);
            ylabel(['y offset = ' num2str(oy, 3)]);
        end

    end
end

%% ------------------------------------------------------------------------
%Functions for random smapling of offset locations - should all return 2
%outputs, offset-x and offset-y, and a single input parameters structure

%Uniform sampling - assume params are |x_max|, |y_max|
function [offset_x, offset_y] = uniform(params)

offset_x = 2*rand*params(1) - params(1);
offset_y = 2*rand*params(2) - params(2);

%Gaussian sampling - assume params are sigma_x, sigma_y
function [offset_x, offset_y] = gaussian(params)

offset_x = randn*params(1);
offset_y = randn*params(2);




  
    
