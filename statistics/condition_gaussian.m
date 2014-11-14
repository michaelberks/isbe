function [c_mean c_covar] = condition_gaussian(mean_in, covar_in, conditions)
%
% Return Value:
%
% Conditioning of a Gaussian is a common statistical technique and is described in most texts on multivariate statistics.
% Here we use the nomenclature of "Multivariate Statistical Methods", by Morrison, 2nd Ed.,
% ISBN: 0-07-043186-8, page 91-93, which is a bit old, but available in the Manchester John Rylands library. The code should
% be transparent enough to see where we are coming from by referring to any good text on multivariate statistics.




% check if we do actually know anything about the model
if isempty(conditions)
% If we get here, we can't condition because we don't have any information about the model.
% We just need to ensure the model is in the natural data space and then return it
% model is already in the natural data space
    c_mean = mean_in;
    c_covar = covar_in;
    return;
end

% first get a list of known dims and unknown dims
known_dims = find(~isnan(conditions));
unknown_dims = find(isnan(conditions));

% condition the components

% Compute the conditioned component:
% get a vector of the unknown dims of the mean for this component
m1 = mean_in(unknown_dims);

% compute C11, the partition of the cov mat for the unknown dims
C11 = covar_in(unknown_dims,unknown_dims);

% compute C12, the partition of the cov mat for the unknown and known dims
C12 = covar_in(unknown_dims,known_dims);

% compute C22, the partition of the cov mat for the known dims
C22 = covar_in(known_dims,known_dims);

% compute the distance of x2 from the known dims -- calculate this only once and reuse in the rest of the loop
x_minus_mu = conditions(known_dims) - mean_in(known_dims);
    
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

    x_solved = linfactor (F, x_minus_mu');
    c_mean = (m1' + (C12 * x_solved))';
    c_covar = C11 - C12 * linfactor (F, C12');
    clear F;
    warning('on', 'MATLAB:nearlySingularMatrix');
else
    clear F;
    pinv_C22 = pinv(C22);
    c_mean = (m1' + (C12 * (pinv_C22 * (x_minus_mu'))))';
    c_covar = C11 - C12 * (pinv_C22 *(C12'));
end

%Make sure c_covar is actually symmetric (it might not be due to numerical
%rounding errors)
c_covar = (c_covar + c_covar')/2;

end