function [conditioned_model] = mb_gmm_condition(aModel, aConditions, aMethod, aPCA_GMMConditioningFuncData)
%
% MB_GMM_CONDITION Condition a Gaussian Mixture Model
%
% This function conditions a GMM, imposing conditions from the natural
% data space. The function copes with regular GMMs, as well as GMMs
% built in a Principal Components space (see note at start of code!!!).
% If conditioning a GMM built in a PC space, the conditions must come from
% the natural data space, and the conditional distribution will be in
% the natural data space (not the full space, as conditions will have
% been applied). If you are working with GMMs in a PC space, then
% ensure you read the notes in this file about providing code to
% condition the components.
%
% Mandatory Arguments:
%
% aModel
%   - The model to condition. This MB_GMM_CONDITION function is written to be able to
%   accept models created by either the MB_K_MEANS_CLUSTERING or MB_CLUSTER_LARGE_DATA_SET
%   (using the MB_K_MEANS_CLUSTERING as the final pass clustering algorithm). However, a Gaussian
%   Mixture Model, generated using any suitable algorithm, can be used here provided that it meets
%   the following specifications: It must be a struct with (at least) the following fields:
%      + 'NumClusters': the number of components (i.e. Gaussians) in the model.
%      + 'ClusterProbs': a vector where ClusterProbs(i) is computed from the number of data points
%			assigned to cluster i, and represents the probability of each cluster.
%      + 'Means': a matrix where Means(i,:) is the mean vector for cluster i.
%      + 'CovMats': a cell array where CovMats{i} returns the covariance matrix for cluster i
%   If the model is built in a PC space, the model struct must contain a field called
%   'pca_info_struct', as returned by st_incremental_pca.
%
% aConditions
%   - This is a vector which has **the same length as the model's dimensionality**, and specifies how
%   the model should be conditioned. The model is conditioned using the values in the vector. If
%   no information is available for a set of dimensions, then this should be specified by setting the
%   corresponding elements in the vector to NaN. So, for example, conditioning a 5-dimensional
%   model with the vector [1.2 NaN 3 3.2 NaN] would result in a 2-dimensional model being returned,
%   representing dimensions 2 and 5, and the model would be conditioned on elements 1, 3 and 4 of
%   the above vector.
%
% Optional Arguments:
%
% 'PCA_GMMConditioningFunc':
%
% aMethod
%   - the method used to condition the model. Can use CJR's original method (default option), which works very well,
%		but always selects the most likely cluster. Or, can use Paul Baggenstoss' method which computes the conditional
%		component probabilities. Note that Baggenstoss' method does not (yet) 	support conditioning GMMs in a PC
%		space -- use the CJR method for this, or edit the code surrounding the Baggenstoss method). NOTE: Having
%		tested both methods (see TEST_GMM_CONDITION), CJR's method seems to produce better results than
%		Baggenstoss' method. Although Baggenstoss' method would seem to be better, computing probabilities
%		for every component, CJR's method produces better results.
%	
%
% Return Value:
%
% MB_GMM_CONDITION returns a new model, which represents the conditioning of the model
% specified by the 'Model' argument (see above). The new model conforms to
% specification detailed above.

% cope with optional arguments.
if nargin < 3
    aMethod = 'CJR';
    aPCA_GMMConditioningFuncData = [];
end

% Conditioning of a Gaussian is a common statistical technique and is described in most texts on multivariate statistics.
% Here we use the nomenclature of "Multivariate Statistical Methods", by Morrison, 2nd Ed.,
% ISBN: 0-07-043186-8, page 91-93, which is a bit old, but available in the Manchester John Rylands library. The code should
% be transparent enough to see where we are coming from by referring to any good text on multivariate statistics.

% work out which method to use
if strcmp(upper(aMethod), 'CJR')
    % check to see if the model is in a PCA space
    if isfield(aModel, 'pca_info_struct')
        % the field exists, but is the model in the PCA space or the natural data space?
        % You're in the wrong place fella, try CJR's function
        [conditioned_model] = ...
            mb_gmm_condition(aModel, aConditions, aMethod, aPCA_GMMConditioningFuncData);
        return;
    else
    %model is in natural space - we can MB's modified verison of CJR's
    %original function
        
        % check if we do actually know anything about the model
        if isempty(aConditions)
        % If we get here, we can't condition because we don't have any information about the model.
        % We just need to ensure the model is in the natural data space and then return it
        % model is already in the natural data space
            conditioned_model = aModel;
            return;
        end
        num_dims_natural_space = length(aModel.Means(1,:));
        % check that the vector of conditions is the same length as the
        % number of dims in the model
        if length(aConditions) ~= num_dims_natural_space
            error('The vector of conditions must have the same dimensionality as the model')
        end
        
        % first get a list of known dims and unknown dims
        known_dims = find(~isnan(aConditions));
        unknown_dims = find(isnan(aConditions));
        
        % get x2, the vector of the known measurements
        x2 = reshape(aConditions(known_dims), 1, []);

        % Set the NumClusters and ClusterProbs fields of conditioned_model so we can use the cjr_gmm_marginalise function on it, 
        % and normalise the ClusterProbs so that they sum to unity (i.e. define a pdf)
        conditioned_model.NumClusters = aModel.NumClusters;
       
        % condition the components
        mahals = zeros(1, aModel.NumClusters); % pre-allocate
%         r_good = Inf;
%         r_bad = 0;
        for i = 1 : aModel.NumClusters
            % First compute the conditioned component:
            % get a vector of the unknown dims of the mean for this component
            m1 = aModel.Means(i, unknown_dims);

            % compute C11, the partition of the cov mat for the unknown dims
            C11 = aModel.CovMats{i}(unknown_dims,unknown_dims);

            % compute C12, the partition of the cov mat for the unknown and known dims
            C12 = aModel.CovMats{i}(unknown_dims,known_dims);

            % compute C22, the partition of the cov mat for the known dims
            C22 = aModel.CovMats{i}(known_dims,known_dims);
            
            % compute the distance of x2 from the known dims -- calculate this only once and reuse in the rest of the loop
            x2_minus_mu = x2 - aModel.Means(i,known_dims);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % The next lines of code take up 98% of processing time. Original method used
            % Moore-Penrose pseudo-inverses however there is much faster method
            % (x10)using linsolve: i.e. factor equations using LU or chol, then solve
            % this needs testing though. Note that cjr's pinv_cov is slower than
            % Matlab's pinv now, maybe due to improvement in inbuilt matlab
            % implementation

            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Old code
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %             
            %            pinv_C22 = pinv(C22);
            %            conditioned_model.Means(i, :) = (m1' + (C12 * (pinv_C22 * (x2_minus_mu'))))';
            %            conditioned_model.CovMats{i} = C11 - C12 * (pinv_C22 *(C12'));
            %            mahals(i) = x2_minus_mu * pinv_C22 * x2_minus_mu';
            %
            %             % MAB: Use \ not linfactor
            %             conditioned_model.Means(i, :) = (m1' + (C12 * (C22 \ ((x2 - m2)'))))';
            %             conditioned_model.CovMats{i} = C11 - C12 * (C22 \ (C12'));
            %             warning('on', 'MATLAB:nearlySingularMatrix');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % New stuff
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
%             F = linfactor(C22);
%             if F.code == 2
%                 rcnd = (min (abs(diag(F.L))) / max(abs(diag(F.L)))) ^ 2;
%             else
%                 rcnd = min (abs (diag (F.U))) / max (abs (diag (F.U))) ;
%             end
%             if rcond(C22) > 1e-9
%                 warning('off', 'MATLAB:nearlySingularMatrix');
%                 % compute the conditional distribution for a non-PCA model
%                 
%                 x_solved = linfactor (F, x2_minus_mu');
%                 conditioned_model.Means(i, :) = (m1' + (C12 * x_solved))';
%                 conditioned_model.CovMats{i} = C11 - C12 * linfactor (F, C12');
%                 mahals(i) = x2_minus_mu* x_solved;
%                 %clear F;
%                 if rcnd < 1e-5
%                     pinv_C22 = pinv(C22);
%                     display(['error, rcnd limit too high, rcnd = ', num2str(rcnd)]);
%                     display(['error in mahals = ', num2str(mahals(i) -  x2_minus_mu * pinv_C22 * x2_minus_mu')]);
%                     display(['correct mahals = ', num2str(x2_minus_mu * pinv_C22 * x2_minus_mu')]);
%                 end
%                 r_good = min(r_good, rcnd);
%                 warning('on', 'MATLAB:nearlySingularMatrix');
%             else
%                 pinv_C22 = pinv(C22);
%                 conditioned_model.Means(i, :) = (m1' + (C12 * (pinv_C22 * (x2_minus_mu'))))';
%                 conditioned_model.CovMats{i} = C11 - C12 * (pinv_C22 *(C12'));
%                 mahals(i) = x2_minus_mu * pinv_C22 * x2_minus_mu';
%                 if rcnd >= 1e-5
%                     x_solved = linfactor (F, x2_minus_mu');
%                     display(['error, rcnd limit too low, rcnd = ', num2str(rcnd)]);
%                     display(['error in mahals = ', num2str(mahals(i) - x2_minus_mu* x_solved)]);
%                     display(['correct mahals = ', num2str(mahals(i))]);
%                 end
%                 r_bad = max(r_bad, rcnd);
%                 %display(['Dist = ', num2str(mahals(i))]);
%             end
            
            try
                F.L = chol(C22, 'lower');
                F.code = 2;
                % if chol works then test condition
                rcnd = (min (abs(diag(F.L))) / max(abs(diag(F.L)))) ^ 2 > 1e-5;
            catch
                rcnd = false; %if chol fails we know matrix is nearly singular so no need to compute condition
            end
            
            if rcnd
                warning('off', 'MATLAB:nearlySingularMatrix');
                % compute the conditional distribution for a non-PCA model
                
                x_solved = linfactor (F, x2_minus_mu');
                conditioned_model.Means(i, :) = (m1' + (C12 * x_solved))';
                conditioned_model.CovMats{i} = C11 - C12 * linfactor (F, C12');
                mahals(i) = x2_minus_mu* x_solved;
                clear F;
%                 if rcnd < 1e-5
%                     pinv_C22 = pinv(C22);
%                     display(['error, rcnd limit too high, rcnd = ', num2str(rcnd)]);
%                     display(['error in mahals = ', num2str(mahals(i) -  x2_minus_mu * pinv_C22 * x2_minus_mu')]);
%                     display(['correct mahals = ', num2str(x2_minus_mu * pinv_C22 * x2_minus_mu')]);
%                 end
%                 r_good = min(r_good, rcnd);
                warning('on', 'MATLAB:nearlySingularMatrix');
            else
                clear F;
                pinv_C22 = pinv(C22);
                conditioned_model.Means(i, :) = (m1' + (C12 * (pinv_C22 * (x2_minus_mu'))))';
                conditioned_model.CovMats{i} = C11 - C12 * (pinv_C22 *(C12'));
                mahals(i) = x2_minus_mu * pinv_C22 * x2_minus_mu';
%                 if rcnd >= 1e-5
%                     x_solved = linfactor (F, x2_minus_mu');
%                     display(['error, rcnd limit too low, rcnd = ', num2str(rcnd)]);
%                     display(['error in mahals = ', num2str(mahals(i) - x2_minus_mu* x_solved)]);
%                     display(['correct mahals = ', num2str(mahals(i))]);
%                 end
%                 r_bad = max(r_bad, rcnd);
                %display(['Dist = ', num2str(mahals(i))]);
            end
            
        
           
        end
%         display(['r_good = ', num2str(r_good), ' r_bad = ', num2str(r_bad)]);
        
        % Now compute the new mixing proportions: the new mixing proportions are
        % proportional to the product of the heights of the clusters in the marginalised
        % distribution just computed (evaluated at x2) and the mixing
        % proportions of the original distribution
        % compute the mahalanobis distance for this component           
        % compute the conditional component probabilities
        
        %MB: divide Mahals by their sum - the absolute value is not
        %important just the proportions and e^-km = (e^-m).(e^k). This will
        %avoid errors due to roundup (i.e. all cluster probs = 0)
        mahals = mahals / sum(mahals);
        conditioned_model.ClusterProbs = aModel.ClusterProbs .* exp(-mahals / 2);
        
%         old_clusterprobs = conditioned_model.ClusterProbs;
%         warning off MATLAB:divideByZero; % may get a divide by zero below, and
%         I don't want to waste time printing warning messages to the screen
        conditioned_model.ClusterProbs = conditioned_model.ClusterProbs ./...
            sum(conditioned_model.ClusterProbs);
%         warning on  MATLAB:divideByZero; % OK, now just get on with things
% 
%         if any(isnan(conditioned_model.ClusterProbs)) || ~fp_compare(sum(conditioned_model.ClusterProbs), 1)
%             % We have computed the component probabilities wrongly
%             % choose a cluster based on its Mahalanobis distance (choose closest)
%             % -- set its prob to 1 and all other probs to zero
%             [dummy smallest_Mahal_index] = min(mahals);
%             conditioned_model.ClusterProbs(:) = 0;
%             conditioned_model.ClusterProbs(smallest_Mahal_index) = 1;
%         end
    end
    
elseif strcmp(upper(aMethod), 'BAGGENSTOSS')
    % convert to the form used by Baggenstoss
    for i = 1 : aModel.NumClusters
        param.modes(i).weight = aModel.ClusterProbs(i);
        param.modes(i).mean = aModel.Means(i, :)';
        R = chol(aModel.CovMats{i});
        param.modes(i).cholesky_covar = R;
        param.features(i).name = i;
        param.features(i).min_std = 0.1;
    end
    
    % call the Baggenstoff code
    param = gmix_condx(param, find(~isnan(aConditions)),...
        find(isnan(aConditions)), aConditions(~isnan(aConditions)));
    
    % convert back to the form used by CJR
    for i = 1 : aModel.NumClusters
        conditioned_model.ClusterProbs(i) = param.modes(i).weight;
        conditioned_model.Means(i, :) = param.modes(i).mean';
        
        % the covariance matrices are returend as their Cholesky decomposition, so
        % convert to the explicit covariance matrices
        conditioned_model.CovMats{i} = (param.modes(i).cholesky_covar' * param.modes(i).cholesky_covar');
    end
else
    error(['Method ' aMethod ' is not supported']);
end