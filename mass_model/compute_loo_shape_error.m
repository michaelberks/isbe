% COMPUTE_LOO_SHAPE_ERROR put mass outlines into shape matrix 
%    [loo_errors] = compute_loo_shape_error(shapes)
%
%    inputs:
%       shapes  - matrix of shapes to apply leave-one-out testing to
%       thresh  - threshold for PCA can either by [0,1] specifying the
%           percantage of variance to keep, or a +ve int spcuifying the
%           number of modes to keep
%    outputs:
%       loo_error - vector of loo RMS errors
%
%    notes: MDL values? other outputs?
%
% Copyright: (C) 2006-2008 University of Manchester
% Author: Michael Berks

function [loo_errors variances] = compute_loo_shape_error(shapes, n_modes)

N = size(shapes, 1);    
M = length(n_modes);

loo_errors = zeros(N, M);
variances = zeros(N, M);

for ii = 1:N
    %remove i-th shape from data and set as unseen test shape
    loo_shapes = shapes([1:ii-1, ii+1:end], :);
    shape_t = shapes(ii, :);
    
    for jj = 1:M
        % Compute model parameters for shapes
        [mean_loo, P_loo, dummy, L_loo] = pca(loo_shapes, n_modes(jj));
        B_loo = P_loo' * (shape_t - mean_loo)';
        
        shape_rg = mean_loo + (P_loo*B_loo)';

        % Compute loo errors a RMS difference between regenerated shape and
        % test shape
        loo_errors(ii,jj) = sqrt(mean(sum((reshape(shape_rg, [], 2) -...
                            reshape(shape_t, [], 2)).^2, 2)));
        variances(ii,jj) = sum(L_loo);
    end
    
end 
    
    
    