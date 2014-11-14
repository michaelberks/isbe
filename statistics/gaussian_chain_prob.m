function [log_p] = gaussian_chain_prob(contour_xy, widths, full_mu, full_sigma, P_known, z_mu, z_sigma, parts_length)
%GAUSSIAN_CHAIN_PROB *Insert a one line summary here*
%   [p] = gaussian_chain_prob(contour_xy, c_mu, c_sigma, parts_length)
%
% Inputs:
%      contour_xy - *Insert description of input variable here*
%
%      c_mu - *Insert description of input variable here*
%
%      c_sigma - *Insert description of input variable here*
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


num_pts = size(contour_xy,1);
num_modes = size(P_known,2);
known_dims = 1:num_modes;
known_dims_z = 1:3*(parts_length-1);
con_dims_z = 3*(parts_length-1)+(1:3);

%Reshape mean segment from z_mu to act as target shape for alignment
mean_segment = reshape(z_mu(1:2*(parts_length-1)), parts_length-1, 2);


for i_pt = parts_length:num_pts
    
    %Get segment and align it to the mean
    segment = contour_xy(i_pt-parts_length+1:i_pt-1,:);
    [segment transform] = align_to_mean(segment, mean_segment);
    segment_widths = widths(i_pt-parts_length+1:i_pt-1);
    z_segment = ([segment(:)' segment_widths(:)'] - z_mu(known_dims_z)) ./ z_sigma(known_dims_z);
    
    %Transform aligned segment into the principal space
    b_segment = z_segment*P_known;
    
    %On first iteration, compute prob of intitial segment
    if i_pt == parts_length
        log_p = log( mvnpdf(b_segment, full_mu(known_dims), full_sigma(known_dims,known_dims)) );
    end
    
    %Condition the Gaussian of the next point using the segment
    conditions = NaN(1,num_modes+3);
    conditions(known_dims) = b_segment;
    [c_mu c_sigma] = condition_gaussian(full_mu, full_sigma, conditions);
    
    %Align the next point using the segment's transform
    next_pt = contour_xy(i_pt,:);
    aligned_pt = transform_pt(next_pt, transform);
    c_data = ([aligned_pt widths(i_pt)] - z_mu(con_dims_z)) ./ z_sigma(con_dims_z);
    
    %Compute log prob of the aligned point under the conditioned Gaussian
    log_p_i = log( mvnpdf(c_data, c_mu, c_sigma) );
    log_p = log_p + log_p_i;    
    
end

function [aligned_segment transform] = align_to_mean(segment, mean_segment)

[~, aligned_segment, transform] = mb_procrustes(mean_segment, segment);

function [aligned_pt] = transform_pt(pt, transform)

aligned_pt = transform.b * pt * transform.T + transform.c(1,:);







