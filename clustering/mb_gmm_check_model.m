function [] = mb_gmm_check_model(model)
%
% Check a GMM to make sure it is valid.
%
% Call this function after bulding a GMM and before
% using it, to make sure that the model is OK.
%
% This function assumes a model of the form built by
% MB_GMM_EM.
%
% This function checks that:
%		+ the covariance matrices are symmetric
%		+ the covariance matrices do not have any negative diagonal elements
%		(i.e. so we are sure they are valid covariances and can use PINV_COV
%		without worrying).
%		+ the mixing proportions sum to unity
%		+ the covariance matrices are all the same size
%		+ the means and covariance matrices are the same size
%		+ the reported number of components is correct
%		(i.e. the same as the number of means and covariance matrices)

% Check to ensure that the covariance matrices are symmetric
% and have non-negative diagonals
for i = 1 : model.NumClusters
	if ~mb_issymmetric(model.CovMats{i})
		figure;
		imagesc(model.CovMats{i});
		error(['Covariance matrix ' num2str(i) ' is not symmetric']);
	end
	
	if any(diag(model.CovMats{i}) < 0)
		error(['Covariance matrix ' num2str(i) ' has one or more negative diagonal elements']);
	end
end

% Check that the mixing proportions sum to unity
if ~fp_compare(sum(model.ClusterProbs), 1)
	error('The mixing proportions do not sum to unity');
end

% Check to ensure that the covariance matrices are all the same size
cov_size = size(model.CovMats{1}, 1); % can assume they are square, since at this point they are symmetric
for i = 2 : model.NumClusters
	if cov_size ~= size(model.CovMats{i}, 1)
		error('The covariance matrices are not all the same size');
	end
end

% Check to ensure that the means and covariance matrices are the same size
if cov_size ~= size(model.Means, 2)
	error('The means and the coavariance matrices differ in size');
end

