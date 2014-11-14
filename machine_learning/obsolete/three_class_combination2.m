function [combination,accuracy] = three_class_combination2(num_clusts,cluster_labels,true_labels)
%THREE_CLASS_COMBINATION2 Given true labels for a three class classification
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
    % Workout the accuracy measure given a combination labelling
    %----------------------------------------------------------------------
    %Work out the total number of combinations for this set of classes
    num_combs = 2*(3^(num_clusts-2));
    
    accuracy = 0;
    
    comb2 = [2 2 3 1 3 1];
    comb3 = [1 3 2 3 1 2];

    for ii = 1:num_combs
        
        %Pre-allocate comb_labels and combination
        comb_labels = uint8(zeros(size(true_labels)));
        this_combination = zeros(1, num_clusts);
        %For each cluster work out class (do separately for columns 1 to 3)
        
        for k = 1:num_clusts
            
            %Work out class for this cluster
            if k == 1
                c = rem(ceil(ii / 2), 3);     
            elseif k == 2
                c = comb2(rem(ii,6)+1);
            elseif k == 3
                c = comb3(rem(ii,6)+1);
            else
                c = rem(ceil(ii / (2 * 3^(k-3))), 3);
            end
            if ~c
                c = 3;
            end
            c = uint8(c);
            %Assign class label
            comb_labels(cluster_labels == uint8(k)) = c;
            this_combination(k) = c;    
        end
        display(num2str(this_combination));
        
        %do accuracy measure of comb_labels vs true_labels
        comb_accuracy = sum(comb_labels == true_labels);
        %display(num2str(comb_accuracy));
        
        %If we've improved accuracy save this combination
        if comb_accuracy > accuracy 
            accuracy = comb_accuracy;
            combination = this_combination;
        end
    end
    
    %Get accuracy as a % - do this now to avoid a division on every combination
    accuracy = accuracy / length(true_labels);
