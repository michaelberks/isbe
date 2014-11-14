function [dist_ratio, dist_1_1, dist_1_2] = knn_2_sample_distance(samp1, samp2, k)
%kNN_2_SAMPLE_DISTANCE for each point in 1 sample, compute the ratio of the
%mean distance to its k-nearest neighbours in its own sample and in the
%second sample
%   [dist_ratio, dist_1_1, dist_1_2] = knn_2_sample_distance(samp1, samp2, k)
%
% Inputs:
%      samp1 - N*d sample of points
%
%      samp2 - M*d sample of points
%
%      k - number of neighbours to use
%
%
% Outputs:
%      dist_ratio -  ratio of within sample 1 distances to distance between
%      classes
%
%      dist_1_1 - within sample 1 distances
%
%      dist_1_2 - distances from sample 1 to sample 2
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 24-Jun-2010
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester 

%Pre-allocate storage for k nearest neighbours
num_samp1 = size(samp1, 1);
d_sorted_samp1 = zeros(num_samp1, k);
d_index_samp1 = zeros(num_samp1, k);

d_sorted_samp2 = zeros(num_samp1, k);
d_index_samp2 = zeros(num_samp1, k);

%Look through points in sample 1
for i = 1:num_samp1
    
    %Compute distances to each other point in sample 1
    Dk = sum(bsxfun(@minus, samp1([1:i-1 i+1:end],:), samp1(i,:)).^2, 2);
    
    %Select the k points with the smallest distances
    if k>1
        [sorted,index] = sort(Dk);
        d_sorted_samp1(i,:) = sorted(1:k);
        d_index_samp1(i,:) = index(1:k);
    else
        [d_sorted_samp1(i,:),d_index_samp1(i,:)] = min(Dk);
    end
    
    %Compute distance to each point in sample 2
    Dk = sum(bsxfun(@minus, samp2, samp1(i,:)).^2, 2);
    
    %Select the k points with the smallest distances
    if k>1
        [sorted,index] = sort(Dk);
        d_sorted_samp2(i,:) = sorted(1:k);
        d_index_samp2(i,:) = index(1:k);
    else
        [d_sorted_samp2(i,:),d_index_samp2(i,:)] = min(Dk);
    end
end

%Compute the mean of the k distances for each point
dist_1_1 = mean(d_sorted_samp1, 2);
dist_1_2 = mean(d_sorted_samp2, 2);

%Compute the ratio of the inter/intra sample distances
dist_ratio = dist_1_2 ./ dist_1_1;
