function [sample, clusters] = mb_gmm_sample(varargin)
%
% MB_GMM_SAMPLE Sample from a Gaussian Mixture Model
%
% MB_GMM_SAMPLE uses the U_PACKARGS interface function
% and so arguments to the function may be passed as name-value pairs
% or as a struct with fields with corresponding names. The field names
% are defined as:
%
% Mandatory Arguments:
%
% 'Model'
%   - The model to sample from. This MB_GMM_SAMPLE function is written to be able to
%   accept models created by either the MB_K_MEANS_CLUSTERING or MB_CLUSTER_LARGE_DATA_SET
%   (using the MB_K_MEANS_CLUSTERING as the final pass clustering algorithm). However, a Gaussian
%   Mixture Model, generated using any suitable algorithm, can be used here provided that it meets
%   the following specifications: It must be a struct with (at least) the following fields:
%      + 'NumClusters': the number of components (i.e. Gaussians) in the model.
%      + 'ClusterProbs': a vector where ClusterProbs(i) is computed from the number of data points
%			assigned to cluster i, and represents the probability of each cluster.
%      + 'Means': a matrix where Means(i,:) is the mean vector for cluster i.
%      + 'CovMats': a cell array where CovMats{i} returns the covariance matrix for cluster i
%
% Optional Arguments:
%
% 'NumSamples'
%   - The number of samples to draw from the distribution described by the model. If not specified,
%   only one sample will be returned.
%
% 'ForceCluster'
%   - Force a particular cluster to be chosen -- i.e. a point is drawn from the specified cluster. 
%   This lets you test specific clusters, e.g. unlikely clusters. Must be an integer in the range
%   1 to Model.NumClusters. Defaults to the empty matrix (i.e. clusters are chosen proportional
%   to their probability).
%
% 'Debug'
%   - Set to true to enable debug functionality. See the code for details of exactly what is
%   available.
%
% Return Value:
%
% MB_GMM_SAMPLE returns sample, which is a matrix containing samples from the distribution described
% by the Gaussian Mixture Model supplied. Each row of the matrix is a sample. The matrix will be a row vector if
% only one sample is requested (the default behaviour).

% first unpack the arguments
args = u_packargs(varargin,... % the user's input
			 'strict', ... % strict mode
			 {'Model'}, ... % the mandatory arguments
			 'NumSamples', 1, ... % the optional arguments
             'ForceCluster', [], ...
             'Debug', (1==0));

% work out the dimensionality of the model, for pre-allocation purposes
sample = zeros(args.NumSamples, length(args.Model.Means(1,:)));
clusters = zeros(args.NumSamples, 1);

% see if the user has force a particular cluster to be chosen
if ~isempty(args.ForceCluster)
    
    for i = 1 : args.NumSamples

        chosen_cluster = args.ForceCluster;
        clusters(i) = chosen_cluster;
        
        % now sample from that cluster
        sample(i,:) = sample_from_normal(args.Model.Means(chosen_cluster,:),...
            args.Model.CovMats{chosen_cluster}, 1);
    end
else
% if not then choose cluster
    
    for i = 1 : args.NumSamples

        %choose a cluster based upon num points assigned to clusters (or whatever measure was used)
        chosen_cluster = mb_sample_based_on_probs('P', args.Model.ClusterProbs, 'Strict', 'notstrict');
        clusters(i) = chosen_cluster;
        
        % now sample from that cluster
        sample(i,:) = sample_from_normal(args.Model.Means(chosen_cluster,:),...
            args.Model.CovMats{chosen_cluster}, 1);
    end
end
% do any debug stuff
if args.Debug
    disp(['Chosen Cluster is: ' num2str(chosen_cluster)]);
end