function [training_data training_labels] = generate_nailfold_vessel_data(varargin)
%SAMPLE_VESSEL_DT_DATA *Insert a one line summary here*
%   [] = sample_mammo_training_data(varargin)
%
% SAMPLE_VESSEL_DT_DATA uses the U_PACKARGS interface function
% and so arguments to the function may be passed as name-value pairs
% or as a struct with fields with corresponding names. The field names
% are defined as:
%
% Mandatory Arguments:
%
% Optional Arguments:
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 17-Aug-2010
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester
warning('off', 'load_uint8:missing_variables');
args = u_packargs(varargin,... % the user's input
    '0', ... % non-strict mode
    {'num_samples', ...
    'image_dir',...
    'fov_mask_dir',...
    'vessel_mask_dir'},... % the mandatory arguments
    'prediction_type', 'detection',...
    'ori_dir', [],...
    'selected_images', [], ...
    'win_size', 3,...
    'num_levels', 4,...
    'feature_type', 'conj',...
    'feature_shape', 'rect',...
    'rotate', 0,...
    'do_max', 0,...
    'decomp_type', 'dt',...
    'bg_ratio', 1,...
    'num_angles', 8,...
    'min_wavelength', 4,...
    'onf', 0.65,...
    'sigma_range', [1 2 4 8],...
    'pca', [],...
    'use_nag', 0);
clear varargin;

num_samples = args.num_samples; %total number of feature vectors to sample

%Copy arguments used in sampling funtion into new structure
switch args.decomp_type
    case 'dt'
        if length(args.num_levels) == 1
            sampling_args.levels = 1:args.num_levels;
        else
            sampling_args.levels = args.num_levels;
        end
        sampling_args.feature_shape = args.feature_shape;
        sampling_args.feature_type = args.feature_type;
        sampling_args.do_max = args.do_max;
        sampling_args.rotate = args.rotate;
        sampling_args.win_size = args.win_size;
    case 'mono'
        sampling_args.win_size = args.win_size;
    case {'g2d', 'g2di', 'clover', 'haar', 'g1d', 'g12d'}
        sampling_args.win_size = args.win_size;
    case 'linop'
        sampling_args.win_size = args.win_size;
        sampling_args.num_levels = args.num_levels;
        sampling_args.num_angles = args.num_angles;
        sampling_args.do_max = args.do_max;
        sampling_args.rotate = args.rotate;
    case 'pixel'
        sampling_args.win_size = args.win_size;
        
    case 'dtg2'
        if length(args.num_levels) == 1
            sampling_args.levels = 1:args.num_levels;
        else
            sampling_args.levels = args.num_levels;
        end
        sampling_args.feature_shape = args.feature_shape;
        sampling_args.feature_type = args.feature_type;
        sampling_args.do_max = args.do_max;
        sampling_args.rotate = args.rotate;
        sampling_args.win_size = args.win_size;
    otherwise
        warning(['decomposition type: ', args.decomp_type, ' not recognised, using DT-CWT']); %#ok
        sampling_args.feature_shape = args.feature_shape;
        sampling_args.feature_type = args.feature_type;
        sampling_args.do_max = args.do_max;
        sampling_args.rotate = args.rotate;
        sampling_args.win_size = args.win_size;
end

%work out size of feature vectors

%See if principal components are specified, and if so, load them. This
%defines the sample feature length
if ~isempty(args.pca);
    %pca should have fields mean and modes
    size_sample_vector = size(args.pca.modes,2);
    sampling_args.pca = args.pca; args = rmfield(args, 'pca');
    
else %Workout feature length from sampling arguments
    
    switch args.decomp_type
        case 'dt'
            size_sample_vector = length(sampling_args.levels)*compute_dt_feature_size(sampling_args);
            
        case 'mono'
            size_sample_vector = 3*args.num_levels*args.win_size^2;
            
        case {'g2d','clover'}
            size_sample_vector = 3*length(args.sigma_range)*args.win_size^2;
            
        case {'haar', 'g1d'}
            size_sample_vector = 2*length(args.sigma_range)*args.win_size^2;
            
        case 'g12d'
            size_sample_vector = 5*length(args.sigma_range)*args.win_size^2;
			
        case 'g2di'
            size_sample_vector = 3*args.sigma_range(2)*args.win_size^2;
            
        case 'linop'
            size_sample_vector = ...
                (args.num_angles - args.do_max*(args.num_angles-1))*args.num_levels*args.win_size^2;
            
        case 'pixel'
            size_sample_vector = args.win_size^2;
            
        case 'dtg2'
            size_sample_vector = ...
                length(sampling_args.levels)*compute_dt_feature_size(sampling_args) +...
                3*length(args.sigma_range)*args.win_size^2;
            
        otherwise
            size_sample_vector = length(sampling_args.levels)*compute_dt_feature_size(sampling_args);
    end    
end

%Create memory storage for data
training_data = zeros(args.num_samples, size_sample_vector);
training_labels = zeros(args.num_samples, 1);

%Get list of all the images, fov masks and vessel masks - we're going to
%assume these all match up for now
image_list = dir([args.image_dir '/*.mat']);
fov_list = dir([args.fov_mask_dir '/*.mat']);
vessel_list = dir([args.vessel_mask_dir '/*.mat']);
switch args.prediction_type
	case {'orientation','linear_regression','logistic_regression','boosted_regression','centre_orientation','mixed_orientation','mixed_centre_orientation'},
		ori_list = dir([args.ori_dir '/*.mat']);
end 

%Check which images are selected - we'll assume if images are selected the
%user has managed to index images within the corrcet range
if isempty(args.selected_images)
    selected_images = 1:length(image_list);
else
    selected_images = args.selected_images;
end
num_images = length(selected_images);

%loop through each image sampling data
curr_sample = 1;
for kk = 1:num_images
    this_im = selected_images(kk);
    
    %Work out the number of samples to take
    num_samples_image = ...
        sample_from_binomial((num_samples + 1 - curr_sample), 1/(num_images+1-kk), 1);
    
    if ~num_samples_image
        continue;
    end
    
    %
    % Get relevant output data for given prediction type
    %----------------------------------------------------------------------
    switch args.prediction_type
        case {'detection','logistic_classification','centre_detection'}
            %Load in vessel mask
            vessel_mask = load_uint8([args.vessel_mask_dir vessel_list(this_im).name]);

            if strcmp(args.prediction_type, 'centre_detection')
                %Create skeletonised vessel mask
                vessel_mask = bwmorph(vessel_mask, 'thin', 'inf');
            end

            %Load in mask of fov region
            fov_mask = u_load([args.fov_mask_dir fov_list(this_im).name]);

            %Check we have enough samples in data
            total_v_pts = sum(vessel_mask(:) & fov_mask(:));
            total_b_pts = sum(~vessel_mask(:) & fov_mask(:));
            num_samples_image = min(num_samples_image, total_v_pts*(1 + args.bg_ratio)); %just in case...
            display(['Samples in image ' num2str(kk) ': ', num2str(num_samples_image)]);

            %Workout number of background and vessel pixels required
            num_v_pts = floor(num_samples_image / (1 + args.bg_ratio));
            num_b_pts = ceil(num_samples_image * args.bg_ratio / (1 + args.bg_ratio));

            %Get random sample of vessel pixels
            v_idx = find(vessel_mask & fov_mask);
            r_idx = randperm(total_v_pts);
            v_idx = v_idx(r_idx(1:num_v_pts));

            %Get random sample of background pixels
            b_idx = find(~vessel_mask & fov_mask);
            r_idx = randperm(total_b_pts);
            b_idx = b_idx(r_idx(1:num_b_pts));

            %Combine the vessel/background samples
            image_idx = [v_idx; b_idx];
            
            %Save the sample labels for this image
            training_labels(curr_sample:num_samples_image+curr_sample-1) = ...
                [true(num_v_pts, 1); false(num_b_pts, 1)];
            
        case {'fov_detection', 'fov_centre_detection'}
            %Load in vessel mask
            vessel_mask = u_load([args.vessel_mask_dir vessel_list(this_im).name]);

            if strcmp(args.prediction_type, 'fov_centre_detection')
                %Create skeletonised vessel mask
                vessel_mask = bwmorph(vessel_mask, 'thin', 'inf');
            end

            %Load in mask of fov region
            fov_mask = u_load([args.fov_mask_dir fov_list(this_im).name]);
            
            %create a disk region around the outer border, we will
            %sample 10% of points from this region
            edge_mask = fov_mask & ~imerode(fov_mask, strel('disk', 32));

            %Check we have enough samples in data
            total_v_pts = sum(vessel_mask(:) & fov_mask(:));
            num_samples_image = min(num_samples_image, total_v_pts*(1 + args.bg_ratio)); %just in case...
            display(['Samples in image ' num2str(kk) ': ', num2str(num_samples_image)]);

            %Workout number of background and vessel pixels required
            num_v_pts = floor(num_samples_image / (1 + args.bg_ratio));
            num_b_pts = ceil(num_samples_image * args.bg_ratio / (1 + args.bg_ratio));
            
            %Split this further into edge and main points
            edge_prop = 0.9;
            num_vm_pts = floor(edge_prop*num_v_pts);
            num_ve_pts = ceil((1-edge_prop)*num_v_pts);
            num_bm_pts = floor(edge_prop*num_b_pts);
            num_be_pts = ceil((1-edge_prop)*num_b_pts);

            %Get random sample of main vessel pixels
            v_idx = find(vessel_mask & fov_mask);
            r_idx = randperm(length(v_idx));
            vm_idx = v_idx(r_idx(1:num_vm_pts));
            
            %Get random sample of edge vessel pixels
            v_idx = find(vessel_mask & edge_mask);
            r_idx = randperm(length(v_idx));
            ve_idx = v_idx(r_idx(1:num_ve_pts));

            %Get random sample of main background pixels
            b_idx = find(~vessel_mask & fov_mask);
            r_idx = randperm(length(b_idx));
            bm_idx = b_idx(r_idx(1:num_bm_pts));
            
            %Get random sample of main background pixels
            b_idx = find(~vessel_mask & edge_mask);
            r_idx = randperm(length(b_idx));
            be_idx = b_idx(r_idx(1:num_be_pts));

            %Combine the vessel/background samples
            image_idx = [vm_idx; ve_idx; bm_idx; be_idx];
            
            %Save the sample labels for this image
            training_labels(curr_sample:num_samples_image+curr_sample-1) = ...
                [true(num_v_pts, 1); false(num_b_pts, 1)];
            
        case 'junction_detection'
            %Load in vessel mask
            vessel_mask = u_load([args.vessel_mask_dir vessel_list(this_im).name]);

            v_thick = bwdist(~vessel_mask);
            [rows cols] = size(vessel_mask);
            xx = repmat(1:cols, rows, 1);
            yy = repmat((1:rows)', 1, cols);
    
            %Create skeletonised vessel mask
            gtt = bwmorph(vessel_mask, 'thin', 'inf');
            
            %Get xy coords of points
            [c_y c_x] = find(gtt);
            num_pts = length(c_x);
            
            %Check if each point is a juntion point
            ws = 5;
            ws2 = floor(ws/2);
            offs = -ws2:ws2;
            count_mask = ones(ws); count_mask(2:ws-1,2:ws-1) = 0;
            junction_mask = false(rows, cols);
            
            for ii = 1:num_pts

                patch = gtt(c_y(ii) + offs, c_x(ii) + offs);
                patch = bwselect(patch, ws2+1, ws2+1, 8) & count_mask;
                label = bwlabel(patch, 8);
                
                if max(label(:)) >= 3
                    rad = v_thick(c_y(ii), c_x(ii));
                    circ = (xx - c_x(ii)).^2 + (yy - c_y(ii)).^2 < rad^2;
                    circ = circ & vessel_mask;
                    junction_mask = junction_mask | circ;
                end
            end

            %Check we have enough samples in data
            total_v_pts = sum(junction_mask(:));
            total_b_pts = sum(~junction_mask(:) & vessel_mask(:));
            num_samples_image = min(num_samples_image, total_v_pts*(1 + args.bg_ratio)); %just in case...
            display(['Samples in image ' num2str(kk) ': ', num2str(num_samples_image)]);

            %Workout number of background and vessel pixels required
            num_v_pts = floor(num_samples_image / (1 + args.bg_ratio));
            num_b_pts = ceil(num_samples_image * args.bg_ratio / (1 + args.bg_ratio));

            %Get random sample of vessel pixels
            v_idx = find(junction_mask);
            r_idx = randperm(total_v_pts);
            v_idx = v_idx(r_idx(1:num_v_pts));

            %Get random sample of background pixels
            b_idx = find(~junction_mask & vessel_mask);
            r_idx = randperm(total_b_pts);
            b_idx = b_idx(r_idx(1:num_b_pts));

            %Combine the vessel/background samples
            image_idx = [v_idx; b_idx];
            
            %Save the sample labels for this image
            training_labels(curr_sample:num_samples_image+curr_sample-1) = ...
                [true(num_v_pts, 1); false(num_b_pts, 1)];
        
        case 'junction_centre_detection'
            %Load in vessel mask
            vessel_mask = u_load([args.vessel_mask_dir vessel_list(this_im).name]);
    
            %Create skeletonised vessel mask
            centre_mask = bwmorph(vessel_mask, 'thin', 'inf');
            
            %Get xy coords of points
            [c_y c_x] = find(centre_mask);
            num_pts = length(c_x);
            
            %Check if each point is a juntion point
            ws = 5;
            ws2 = floor(ws/2);
            offs = -ws2:ws2;
            count_mask = ones(ws); count_mask(2:ws-1,2:ws-1) = 0;
            junction_mask = false(size(vessel_mask));
            
            for ii = 1:num_pts

                patch = centre_mask(c_y(ii) + offs, c_x(ii) + offs);
                patch = bwselect(patch, ws2+1, ws2+1, 8) & count_mask;
                label = bwlabel(patch, 8);
                junction_mask(c_y(ii) + offs, c_x(ii) + offs) = max(label(:)) >= 3;
                
            end

            %Check we have enough samples in data
            total_v_pts = sum(junction_mask(:));
            total_b_pts = sum(~junction_mask(:) & vessel_mask(:));
            num_samples_image = min(num_samples_image, total_v_pts*(1 + args.bg_ratio)); %just in case...
            display(['Samples in image ' num2str(kk) ': ', num2str(num_samples_image)]);

            %Workout number of background and vessel pixels required
            num_v_pts = floor(num_samples_image / (1 + args.bg_ratio));
            num_b_pts = ceil(num_samples_image * args.bg_ratio / (1 + args.bg_ratio));

            %Get random sample of vessel pixels
            v_idx = find(junction_mask);
            r_idx = randperm(total_v_pts);
            v_idx = v_idx(r_idx(1:num_v_pts));

            %Get random sample of background pixels
            b_idx = find(~junction_mask & vessel_mask);
            r_idx = randperm(total_b_pts);
            b_idx = b_idx(r_idx(1:num_b_pts));

            %Combine the vessel/background samples
            image_idx = [v_idx; b_idx];
            
            %Save the sample labels for this image
            training_labels(curr_sample:num_samples_image+curr_sample-1) = ...
                [true(num_v_pts, 1); false(num_b_pts, 1)];
                
        case {'orientation','linear_regression','logistic_regression','boosted_regression', 'centre_orientation'}
            %Load in vessel mask
            vessel_mask = u_load([args.vessel_mask_dir vessel_list(this_im).name]);
            
            if strcmp(args.prediction_type, 'centre_orientation')
                %Create skeletonised vessel mask
                vessel_mask = bwmorph(vessel_mask, 'thin', 'inf');
            end

            %Load in orientation_map
            ori_map = load_uint8([args.ori_dir ori_list(this_im).name]);

            %Check we have enough samples in data
            total_v_pts = sum(vessel_mask(:));
            num_samples_image = min(num_samples_image, total_v_pts); %just in case...
            display(['Samples in image ' num2str(kk) ': ', num2str(num_samples_image)]);

            %Get random sample of vessel pixels
            v_idx = find(vessel_mask);
            r_idx = randperm(total_v_pts);
            image_idx = v_idx(r_idx(1:num_samples_image));
            
            %Save the sample orientations for this image
            training_labels(curr_sample:num_samples_image+curr_sample-1) = ...
                ori_map(image_idx);
            
        case {'mixed_orientation', 'mixed_centre_orientation'}
            
            %Load in vessel mask
            vessel_mask = u_load([args.vessel_mask_dir vessel_list(this_im).name]);
            
            if strcmp(args.prediction_type, 'mixed_centre_orientation')
                %Create skeletonised vessel mask
                vessel_mask = bwmorph(vessel_mask, 'thin', 'inf');
            end

            %Load in orientation ground truth
            ori_gt = load([args.ori_dir ori_list(this_im).name]);
            
            %Check we have enough samples in data
            total_v_pts = sum(vessel_mask(:));
            num_samples_image = min(num_samples_image, total_v_pts); %just in case...
            display(['Samples in image ' num2str(kk) ': ', num2str(num_samples_image)]);

            %Get random sample of vessel pixels
            v_idx = find(vessel_mask);
            r_idx = randperm(total_v_pts);
            image_idx = v_idx(r_idx(1:num_samples_image));
            
            %now randomly select one of the mixed orientations for each
            %point
            %DEBUG: figure; imgray(vessel_mask);
            	
            ori_sample = zeros(num_samples_image,1);
            for ii = 1:num_samples_image
                mixed_oris = ori_gt.mixed_oris{ori_gt.mixed_idx(image_idx(ii))};
                ori_sample(ii) = mixed_oris(ceil(rand*length(mixed_oris)));
                
                %DEBUG: [y x] = ind2sub(size(vessel_mask), image_idx(ii));
                %DEBUG: quiver(ones(length(mixed_oris),1)*x,ones(length(mixed_oris),1)*y,real(mixed_oris),-imag(mixed_oris), 'color','r');
                %DEBUG: quiver(ones(length(mixed_oris),1)*x,ones(length(mixed_oris),1)*y,-real(mixed_oris),imag(mixed_oris), 'color','r');
            end
            
            %Save the sample orientations for this image
            training_labels(curr_sample:num_samples_image+curr_sample-1) = ...
                ori_sample;
            
    end   
            
    %Load in retinogram
    nailfold = u_load([args.image_dir image_list(this_im).name]);
            
    %Convert to row/col
    [rows cols] = ind2sub(size(nailfold), image_idx);
     
    %
    % 4. Sample the decomposition coefficients
    %----------------------------------------------------------------------
    
    switch args.decomp_type
        case 'dt' 
            %Compute dual-tree transform of image
            dt = compute_dual_tree(nailfold, sampling_args.levels(end), args.use_nag);

            %Sample DT coefficients from specified rows and cols according to
            %sampling arguments
            training_data(curr_sample:num_samples_image+curr_sample-1,:) = ...
                sample_dt_data(dt, rows, cols, sampling_args);
            clear dt;

        case 'mono'
            [local_amp local_phase local_ori] = monogenic(nailfold, args.num_levels, args.min_wavelength, 2, args.onf, 1);
            clear nailfold;
            training_data(curr_sample:num_samples_image+curr_sample-1,:) = ...
                sample_monogenic_data(local_amp, local_phase, local_ori, rows, cols, sampling_args);
            clear local_amp local_phase;

        case 'g2d'
            g2d_responses = compute_gaussian_2nd_derivatives(nailfold,  args.sigma_range);
            clear nailfold;
            training_data(curr_sample:num_samples_image+curr_sample-1,:) = ...
               sample_g2d_data(...
               g2d_responses(:,:,:,1),... 
               g2d_responses(:,:,:,2),...
               g2d_responses(:,:,:,3), rows, cols, sampling_args);
           clear g2d_responses;
        case 'clover'
            clover_responses = compute_clover_responses(nailfold, args.sigma_range);
            clear nailfold;
            training_data(curr_sample:num_samples_image+curr_sample-1,:) = ...
               sample_g2d_data(...
               clover_responses(:,:,:,1),... 
               clover_responses(:,:,:,2),...
               clover_responses(:,:,:,3), rows, cols, sampling_args);
           clear clover_responses;
        case 'haar'
            haar_responses = compute_haar_responses(nailfold, args.sigma_range);
            clear nailfold;
            training_data(curr_sample:num_samples_image+curr_sample-1,:) = ...
               sample_haar_data(...
               haar_responses(:,:,:,1),... 
               haar_responses(:,:,:,2), rows, cols, sampling_args);
           clear haar_responses;
        case 'g2di'       
            g2d_responses = compute_gaussian_2nd_derivatives_d(nailfold,  args.sigma_range(1), args.sigma_range(2));
            clear nailfold;
            training_data(curr_sample:num_samples_image+curr_sample-1,:) = ...
                sample_g2d_data_d(g2d_responses, rows, cols, sampling_args);
            clear g2d_responses;
        case 'g1d'
            g1d_responses = compute_gaussian_1st_derivatives(nailfold,  args.sigma_range);
            clear nailfold;
            training_data(curr_sample:num_samples_image+curr_sample-1,:) = ...
               sample_g1d_data(...
               g1d_responses(:,:,:,1),... 
               g1d_responses(:,:,:,2), rows, cols, sampling_args);
           clear g1d_responses;

       case 'g12d'
            g1d_responses = compute_gaussian_1st_derivatives(nailfold,  args.sigma_range);
            g2d_responses = compute_gaussian_2nd_derivatives(nailfold,  args.sigma_range);
            clear nailfold;

            training_data(curr_sample:num_samples_image+curr_sample-1,:) = [...
               sample_g1d_data(...
               g1d_responses(:,:,:,1),... 
               g1d_responses(:,:,:,2), rows, cols, sampling_args)... 
               sample_g2d_data(...
               g2d_responses(:,:,:,1),... 
               g2d_responses(:,:,:,2),...
               g2d_responses(:,:,:,3), rows, cols, sampling_args)];

           clear g2d_responses g1d_responses;

        case 'linop'
            training_data(curr_sample:num_samples_image+curr_sample-1,:) = ...
                sample_linop_data(nailfold, rows, cols, sampling_args);
            clear nailfold;

        case 'pixel'
            %Sample pixels from image
            training_data(curr_sample:num_samples_image+curr_sample-1,:) = ...
                sample_pixel_data(double(nailfold), rows, cols, sampling_args);

        case 'dtg2'
            %Compute dual-tree transform of image
            dt = compute_dual_tree(nailfold, sampling_args.levels(end), args.use_nag);
            g2d_responses = compute_gaussian_2nd_derivatives(nailfold,  args.sigma_range);
            clear nailfold;

            %Sample DT coefficients from specified rows and cols according to
            %sampling arguments
            training_data(curr_sample:num_samples_image+curr_sample-1,:) = ...
                [sample_dt_data(dt, rows, cols, sampling_args) ...
                sample_g2d_data(...
                g2d_responses(:,:,:,1),... 
                g2d_responses(:,:,:,2),...
                g2d_responses(:,:,:,3), rows, cols, sampling_args)];
            clear dt g2d_responses;

        otherwise
            %Compute dual-tree transform of image
            dt = compute_dual_tree(nailfold, sampling_args.levels(end), args.use_nag);

            %Sample DT coefficients from specified rows and cols according to
            %sampling arguments
            training_data(curr_sample:num_samples_image+curr_sample-1,:) = ...
                sample_dt_data(dt, rows, cols, sampling_args);
            clear dt;
    end
    
    %Update the current sample count
    curr_sample = curr_sample + num_samples_image;
end


