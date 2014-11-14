function [combination,accuracy] = three_class_combination(num_clusts,cluster_labels,true_labels)
%THREE_CLASS_COMBINATION Given true labels for a three class classification
%problem and cluster assignments for each point, work out the best of
%combination of cluster labels to class labels
%   [combination,accuracy] = three_class_combination(num_clusts,cluster_labels,true_labels)
%
% Inputs:
%      num_clusts- number of clusters
%
%      cluster_labels- cluster label for each data point
%
%      true_labels- true class element for each data point
%
%
% Outputs:
%      combination- k x 1 vector of class labels for each cluster for the
%      best combination
%
%      accuracy- accuracy of the best combination
%
%
% Example:
%
% Notes: This functions works up to approx 12 clusters with a reasonable
% number of data (e.g. 10^5). Above that memory issue may become a problem
% for storing the combinations
%
% See also:
%
% Created: 02-Nov-2009
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester

%Make things uint8 to save memory
    cluster_labels = uint8(cluster_labels(:));
    true_labels = uint8(true_labels(:));
    
    %----------------------------------------------------------------------
    %First make a list of all possible combinations
    %----------------------------------------------------------------------
    %Basic set of combinations for three classes
    comb3 = uint8([ 1 2 3; 1 3 2; 2 1 3; 2 3 1; 3 1 2; 3 2 1]);
    
    %Work out the total number of combinations for this set of classes
    num_combs = 2*(3^(num_clusts-2));
    
    %Pre-allocate combinations
    combinations = zeros(num_combs, num_clusts);
    
    %Repeat the basic 3-way combinations in the first 3 columns
    combinations(:,1:3) = repmat(comb3, 3^(num_clusts-3), 1);

    %Add in the class label column for each additional cluster
    for k = 4:num_clusts
        reps1 = 2*(3^(k-3));
        reps2 = 3^(num_clusts-k);
        combinations(:,k) = repmat(kron([1 2 3]', ones(reps1, 1)), reps2,1);
    end
    
    %----------------------------------------------------------------------
    % Could save these combinations and load in for later?
    %----------------------------------------------------------------------
    
    %----------------------------------------------------------------------
    % Now workout the accuracy measure given a combination labelling
    %----------------------------------------------------------------------
    accuracy = 0;
    
    for ii = 1:num_combs     
        %Get labels for this combination
        comb_labels = combinations(ii,cluster_labels)';
        
        %do accuracy measure of comb_labels vs true_labels
        comb_accuracy = sum(comb_labels == true_labels);
        %display(num2str(comb_accuracy));
        
        %If we've improved accuracy save this combination
        if comb_accuracy > accuracy 
            accuracy = comb_accuracy;
            combination = combinations(ii,:);
        end
    end
    
    %Get accuracy as a % - do this now to avoid a division on every combination
    accuracy = accuracy / length(true_labels);
