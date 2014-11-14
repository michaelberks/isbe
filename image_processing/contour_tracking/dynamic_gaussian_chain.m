function [best_idx, min_E, best_xy, best_w] = dynamic_gaussian_chain(states_xy, states_w, full_mu, full_sigma, P_known, z_mu, z_sigma, parts_length)
%DYNAMIC_GAUSSIAN_CHAIN *Insert a one line summary here*
%   [p] = dynamic_gaussian_chain(contour_xy, full_mu, full_sigma, parts_length)
%
% Inputs:
%      contour_xy - *Insert description of input variable here*
%
%      full_mu - *Insert description of input variable here*
%
%      full_sigma - *Insert description of input variable here*
%
%      parts_length - *Insert description of input variable here*
%
%
% Outputs:
%      p - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 04-Apr-2013
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

num_modes = size(P_known,2);
known_dims = 1:num_modes;
known_dims_z = 1:3*(parts_length-1);
con_dims_z = 3*(parts_length-1)+(1:3);

%Reshape mean segment from z_mu to act as target shape for alignment
mean_segment = reshape(z_mu(1:2*(parts_length-1)), parts_length-1, 2);



%known_mu = full_mu(1:num_modes);
%known_sigma = full_sigma(1:num_modes,1:num_modes);

num_pts = length(states_xy);

E_values = cell(num_pts,1);
E_idx = cell(num_pts,1);

%we could work this out as we go, but easier to run through and comput this
%once
num_states = zeros(1,num_pts);
for i_pt = 1:num_pts
    num_states(i_pt) = size(states_xy{i_pt},1);
end

E_values{3} = zeros(num_states(1:3));

%Compute static probs for all combinations of the first three points
for v1 = 1:num_states(1)
    xy_1 = states_xy{1}(v1,:);
    w1 = states_w{1}(v1);
    
    for v2 = 1:num_states(2)
        xy_2 = states_xy{2}(v2,:);
        w2 = states_w{2}(v2);
        
        for v3 = 1:num_states(3)
            xy_3 = states_xy{3}(v3,:);
            w3 = states_w{3}(v3);
            
            segment = [xy_1; xy_2; xy_3];
            [segment] = align_to_mean(segment, mean_segment);
            z_segment = ([segment(:)' w1 w2 w3] - z_mu(known_dims_z)) ./ z_sigma(known_dims_z); 
            b_segment = z_segment*P_known;
            
            E_values{3}(v1, v2, v3) = log( mvnpdf(b_segment, full_mu(known_dims), full_sigma(known_dims,known_dims)) );
        end
    end
end
                
%--------------------------------------------------------------------------
%Forward pass
conditions = NaN(1,num_modes+3);                    
for i_pt = 4:num_pts
    
    E_values{i_pt} = zeros(num_states(i_pt-[2 1 0]));
    
    for v4 = 1:num_states(i_pt)
        xy_4 = states_xy{i_pt}(v4,:);
        w4 = states_w{i_pt}(v4);
        
        for v3 = 1:num_states(i_pt-1)
            xy_3 = states_xy{i_pt-1}(v3,:);
            w3 = states_w{i_pt-1}(v3);
            
            for v2 = 1:num_states(i_pt-2)
                xy_2 = states_xy{i_pt-2}(v2,:);
                w2 = states_w{i_pt-2}(v2);
                
                shape_log_probs = zeros(num_states(i_pt-3),1);
                
                for v1 = 1:num_states(i_pt-3)
                    xy_1 = states_xy{i_pt-3}(v1,:);
                    w1 = states_w{i_pt-3}(v1);
                    
                    segment = [xy_1; xy_2; xy_3];
                    [segment transform] = align_to_mean(segment, mean_segment);
                    z_segment = ([segment(:)' w1 w2 w3] - z_mu(known_dims_z)) ./ z_sigma(known_dims_z); 
                    b_segment = z_segment*P_known;
                    
                    %Condition the Gaussian of the next point using the segment
                    conditions(known_dims) = b_segment;
                    [c_mu c_sigma] = condition_gaussian(full_mu, full_sigma, conditions);

                    %Align the next point using the segment's transform
                    aligned_pt = transform_pt(xy_4, transform);
                    c_data = ([aligned_pt w4] - z_mu(con_dims_z)) ./ z_sigma(con_dims_z);
    
                    %Compute log prob of the aligned point under the conditioned Gaussian                
                    shape_log_probs(v1) = ...
                        log( mvnpdf(c_data, c_mu, c_sigma) ) + ...
                        E_values{i_pt-1}(v1, v2, v3)';
                end
                
                [min_slp min_idx] = max(shape_log_probs);
                
                E_values{i_pt}(v2, v3, v4) = min_slp;
                E_idx{i_pt}(v2, v3, v4) = min_idx;
            end
        end
    end
end

 %Compute the minimal energy of the final point and the state that achieved it
[min_E, final_idx] = max(E_values{num_pts}(:));

[final_idx1 final_idx2 final_idx3] = ind2sub(num_states(num_pts-[2 1 0]), final_idx);

%------------------------------------------------------------------------
% backward pass
best_idx = zeros(num_pts,1);
best_idx(num_pts) = final_idx3;
best_idx(num_pts-1) = final_idx2;
best_idx(num_pts-2) = final_idx1;

%Working backwards through the points...
for i_pt = num_pts-3:-1:1
    best_idx(i_pt) = E_idx{i_pt+3}(best_idx(i_pt+1), best_idx(i_pt+2), best_idx(i_pt+3));
end

if nargout > 2
    best_xy = zeros(num_pts,2);
    best_w = zeros(num_pts,1);
    for i_pt = 1:num_pts
        best_xy(i_pt,:) = states_xy{i_pt}(best_idx(i_pt),:);
        best_w(i_pt,:) = states_w{i_pt}(best_idx(i_pt));
    end
end

%---------------------------------------------------------------------------------------------
function [aligned_segment transform] = align_to_mean(segment, mean_segment)

[~, aligned_segment, transform] = mb_procrustes(mean_segment, segment);

%---------------------------------------------------------------------------------------------
function [aligned_pt] = transform_pt(pt, transform)

aligned_pt = transform.b * pt * transform.T + transform.c(1,:);

% %---------------------------------------------------------------------------------------------
% function log_p = compute_conditonal_prob(segment, next_pt, mean_segment, P_known, z_mu, z_sigma, full_mu, full_sigma)
% 
%     [b_segment transform] = commpute_b_segment(segment, mean_segment, P_known, z_mu, z_sigma);
% 
%     %Condition the Gaussian of the next point using the segment
%     conditions = [b_segment NaN NaN];
%     [c_mu c_sigma] = condition_gaussian(full_mu, full_sigma, conditions);
%     
%     %Align the next point using the segment's transform
%     aligned_pt = transform_pt(next_pt, transform);
%     
%     %Compute log prob of the aligned point under the conditioned Gaussian
%     log_p = log( mvnpdf(aligned_pt, c_mu, c_sigma) );

% %---------------------------------------------------------------------------------------------    
% function log_p = compute_static_prob(segment, mean_segment, P_known, z_mu, z_sigma, known_mu, known_sigma)
% 
%     b_segment = commpute_b_segment(segment, mean_segment, P_known, z_mu, z_sigma);
%     log_p = log( mvnpdf(b_segment, known_mu, known_sigma) );

% %---------------------------------------------------------------------------------------------
% function [b_segment transform] = commpute_b_segment(segment, mean_segment, P_known, z_mu, z_sigma)
% 
%     %Form segment and align it to the mean
%     [segment transform] = align_to_mean(segment, mean_segment);
%     z_segment = (segment(:)' - z_mu) ./ z_sigma;
%     
%     %Transform aligned segment into the principal space
%     b_segment = z_segment*P_known;
    
    