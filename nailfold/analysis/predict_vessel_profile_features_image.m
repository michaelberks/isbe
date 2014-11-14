function [prediction_map vessel_centre_states vessel_centre_scores] = ...
    predict_vessel_profile_features_image(vessel_patch, vessel_centre, apex_width, rf_predictor, varargin)
%EXTRACT_VESSEL_PROFILE_FEATURES *Insert a one line summary here*
%   [profile_features] = extract_vessel_profile_features(outer_edge, inner_edge, vessel_centre)
%
% Inputs:
%      outer_edge, inner_edge - *Insert description of input variable here*
%
%      vessel_centre - *Insert description of input variable here*
%
%
% Outputs:
%      profile_features - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 17-Apr-2013
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

args = u_packargs(varargin,... % the user's input
    '0', ... % strict mode
    'num_profile_pts', 5,...
    'theta_offsets', pi*[-15 0 15]/180,...
    'scale_offsets', [0.75 1 1.25],...
    'num_displacements', 15,...
    'max_displacement', [],...
    'sigma_vote', 2,...
    'plot', 0,...
    'debug', 0,...
    'quiet', 0);
 clear varargin;

num_pts = size(vessel_centre,1);
num_angles = length(args.theta_offsets);
num_scales = length(args.scale_offsets);

num_samples = args.num_displacements*num_angles;
feature_size = 2*num_scales*args.num_profile_pts;

%Compute edge normals from inner and outer edges
%Project normals out from path, and use profiles to detect outer edge
normal_xy = compute_spline_normals(vessel_centre);

%Compute Gaussian (and its hilbert pair) responses
%Use the half the mean apex width as the central scale
[g2d_responses] = compute_gaussian_2nd_derivatives(vessel_patch, apex_width*args.scale_offsets/2);
[h2d_responses] = compute_hilbert_2nd_derivatives_sep(vessel_patch, apex_width*args.scale_offsets/2);

%Also use the mean vessel width to set the max displacement if not supplied
%by user
if isempty(args.max_displacement)
    args.max_displacement = apex_width;
end

%--------------------------------------------------------------------------
%Leave out the mask stuff for now
%Make an estimated mask of the vessel based on the initial annotation and
%the appex width
% initial_vessel_mask = false(size(vessel_patch));
% vessel_idx = sub2ind(size(vessel_patch), round(vessel_centre(:,2)), round(vessel_centre(:,1)));
% initial_vessel_mask(vessel_idx) = 1;
% initial_vessel_mask = imdilate(initial_vessel_mask, strel('disk', floor(0.8*apex_width2)));
% normal_mask = false(num_pts, args.num_displacements, num_scales, num_angles, args.num_profile_pts);
%--------------------------------------------------------------------------
    
%Set up profile sample positions
profile_pos = linspace(-apex_width/2, apex_width/2, args.num_profile_pts);

%Set up profile sample positions
sample_displacements = linspace(-args.max_displacement, args.max_displacement, args.num_displacements);
    
prediction_map = zeros(size(vessel_patch));
vessel_centre_states = cell(num_pts,1);
vessel_centre_scores = cell(num_pts,1);

if args.plot
    figure; 
    subplot(1,2,1); imgray(vessel_patch);
    plot(vessel_centre(:,1), vessel_centre(:,2));
end

all_dis = [];
all_pred = [];

g = gaussian_filters_1d(args.sigma_vote);
len = length(g);
h_pad = zeros(1, len+1);
v_pad = zeros(len,1);
G = g' * g;
G_00 = [G v_pad; h_pad];
G_10 = [v_pad G; h_pad];
G_01 = [h_pad; G v_pad];
G_11 = [h_pad; v_pad G];
vote_region = (1:len+1) - (len+1)/2;

for i_pt = 1:num_pts
        
    
    %Pre-allocate outputs
    profile_features = zeros(num_samples, feature_size);
    profile_centres = zeros(num_samples, 2);
    profile_vecs = zeros(num_samples, 2);
    profile_dis = zeros(num_samples, 1);
    prediction_map_i = zeros(size(vessel_patch));

    %Compute angle from normal vector
    theta = atan2(-normal_xy(i_pt,2), normal_xy(i_pt,1));
    sample_idx = 1;
    for i_dis = 1:args.num_displacements       
        for theta_offset = args.theta_offsets
            %Store the responses and their associated displacement in the
            %containers
            profile_features(sample_idx,:) = extract_vessel_profile_features(...
                g2d_responses, h2d_responses, vessel_centre(i_pt,:), profile_pos,...
                theta+theta_offset, args.scale_offsets, sample_displacements(i_dis));

            profile_vecs(sample_idx,:) = [cos(theta+theta_offset) -sin(theta+theta_offset)];
            profile_centres(sample_idx,:) = vessel_centre(i_pt,:) +...
                sample_displacements(i_dis)*profile_vecs(sample_idx,:);
            profile_dis(sample_idx,:) = sample_displacements(i_dis);
            
            %Increment the sample index
            sample_idx = sample_idx + 1;
        end        
    end
    predicted_displacements = random_forest_reg_predict(rf_predictor, profile_features, 0);
    predicted_centres = profile_centres - bsxfun(@times, profile_vecs, predicted_displacements);
    
    for i_pr = 1:size(predicted_centres,1)
        predicted_pixel = floor(predicted_centres(i_pr,:));
        subpixel_offset = predicted_centres(i_pr,:) - predicted_pixel;
        x1 = subpixel_offset(1);
        x0 = 1 - x1;
        y1 = subpixel_offset(2);
        y0 = 1 - y1;
        G_i = x0*y0*G_00 + x1*y0*G_10 + x0*y1*G_01 + x1*y1*G_11;
        vote_region_r = predicted_pixel(2)+vote_region;
        vote_region_c = predicted_pixel(1)+vote_region;
        prediction_map_i(vote_region_r, vote_region_c) = ...
            prediction_map_i(vote_region_r, vote_region_c) + G_i;
    end
    [maxima_pos, maxima_vals] = local_image_maxima(prediction_map_i, 10, [], 0);
    maxima_pos(4:end,:) = [];
    maxima_vals(4:end) = [];
    vessel_centre_states{i_pt} = maxima_pos;
    vessel_centre_scores{i_pt} = maxima_vals;
    
    
    prediction_map = prediction_map + prediction_map_i;
    
    %predicted_centres_idx = sub2ind(size(vessel_patch), round(predicted_centres(:,2)), round(predicted_centres(:,1)));
    %prediction_map(predicted_centres_idx) = prediction_map(predicted_centres_idx) + 1;
    
    if args.debug
        figure; imgray(prediction_map_i);
        if size(maxima_pos,1) > 3
            plot(maxima_pos(1:3,1), maxima_pos(1:3,2), 'bx');
        else
            plot(maxima_pos(:,1), maxima_pos(:,2), 'rx');
        end
        
        all_dis = [all_dis; profile_dis];%#ok
        all_pred = [all_pred; predicted_displacements]; %#ok
    end
    if args.plot
        plot(predicted_centres(:,1), predicted_centres(:,2), 'x');
    end
    
        
end

if args.plot   
    subplot(1,2,2); imgray(prediction_map);
    plot(vessel_centre(:,1), vessel_centre(:,2), 'g.');
end


