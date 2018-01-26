function error_struct = ...
    orientation_errors(prediction_errors, centre_errors, n_sets, n_per_set)
% Compute median error over a set (or several bootstraps of a set) of
% errors

if (nargin==0)
    % dummy mode - return empty error structure
    error_struct = struct('all', NaN,...
                          'centre', NaN);
    return
end

if ~exist('n_sets','var'), n_sets = 1; end
if ~exist('n_per_set','var'), n_per_set = inf; end

if (n_sets==1)
    % compute over all measurements
    all_vec = median(abs(prediction_errors));
    centre_vec = median(abs(centre_errors));
else
    % sample n_per_set items n_sets times
    all_vec = zeros(1,n_sets);
    centre_vec = zeros(1,n_sets);

    for i = 1:n_sets
        n_samples = min(n_per_set,length(prediction_errors));
        inds = ceil(rand(1,n_samples)*length(prediction_errors));
        all_vec(i) = median(abs(prediction_errors(inds)));

        n_samples = min(n_per_set,length(centre_errors));
        inds = ceil(rand(1,n_samples)*length(centre_errors));
        centre_vec(i) = median(abs(centre_errors(inds)));
    end
end

error_struct = struct('all', all_vec,...
                      'centre', centre_vec);
                       
                       
