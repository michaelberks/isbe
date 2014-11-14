function val = mab_mdl_objective_function(points, argument_list)
%MAB_MDL_OBJECTIVE_FUNCTION Calculate `goodness' of a parameterisation using MDL
%  VAL = MAB_MDL_OBJECTIVE_FUNCTION(POINTS, {COVARIANCE, range, covariance_integral, weighting, connections})
%  
%  inputs:
%    POINTS - (DxNxE) matrix of reparameterised points.  Sort of
%             optional, if not provided then COVARIANCE must be given,
%             but POINTS is used for calculating the data's range.
% argument_list - ordered list of arguments
%    COVARIANCE - (ExE) Optional.  The covariance matrix.  If not
%                provided then it is computed from POINTS.
%    RANGE - (D) a vector giving the data's range.  This must be given iff POINTS is not.
%  outputs:
%    VAL - The `compactness' of the points
%
%  N - number of points
%  E - number of examples
%  References:
%
%  Rhodri Davies, Learning Shape: Optimal Models for Analysing Natural
%  Variability, PhD thesis, 2002.
%
%  See also M_MDL_OBJECTIVE_FUNCTION.

% u_packargs is too inefficient
% use strict, ordered list
% covariange, range, covariance_integral, weighting, connections, clove

if nargin < 2
    argument_list = {[], [], []};
end
covariance = argument_list{1};
range = argument_list{2};
nModes = argument_list{3};
if isempty(points)
    if isempty(covariance) || isempty(range)
        error('MDL:NoPoints', 'Range and Covariance Matrix must be given if points are not');
    end
    nExamples = size(covariance, 1);
    R = max(range);    
else
	[nExamples, nPoints] = size(points);
    points = points - repmat(mean(points), nExamples, 1);
    K = (points * points') / nPoints;
    covariance = K;
  
    % r is the range of the data in data-space.
    % NOTE: this should really be the original data, but it makes bugger all difference in practice
    r = max(points(:)) - min(points(:));
    
    % R is the maximum possible range of Y in parameter-space.
    R = r * sqrt(nPoints);
end

% RHD: I have a suspicion that R is not quite right.

n_s = nExamples;

if isempty(nModes);
    nModes = nExamples - 1;
end

[v,d] = eig(covariance);
d = diag(d)';

[sorted_d,i] = sort(d);
i = fliplr(i);
i = i(1:(nExamples-1));

% eq 3.6 says to multiply by n_p but how does this make sense for
% the continuous case?  Besides, it means that the values of the
% variance depend on the number of points.
%
% We (conceptually) multiply by n_p when reconstructing instead, so
% variance and Y are not the same as those returned by st_pca, but
% they are the same no matter what n_p is and indeed are the same as
% the continuous case.
variance = d(i);

if any(variance < -eps)
    warning('MDL:NegativeVariance', 'Negative variance')
end
v(:, variance < 0) = 0;
variance(variance< 0) = 0;

% p^m', eq. 3.11
pm_p = v(:, i);

% This is not in Davies 2002, but is in eg. Kotcheff & Taylor 1998
% eq 11 and in st_pca.
div = sqrt(variance * nExamples);
div(div==0) = 1; % Avoid divide by zero warnings
%pm_p = pm_p ./ repmat(div, n_s, 1);
for i=1:n_s,
    pm_p(i,:) = pm_p(i,:) ./ div;
end

% If you do the algebra then this is equivalent to the usual
% formulation.  Y is the parameter matrix of the model's training
% set.
Y = nExamples * covariance * pm_p;

% Why divided by 2?
sigma_max = R/2;

%Delta_min = .01;
%Delta_max = 2;
%tol = .1;
Delta_min = 0.0001;
Delta_max = 0.002;
tol = 0.001;

if sigma_max < Delta_max
    warning('MDL:SmallRange', 'Data has a very small range, Delta will dominate it.');
end

if any(abs(Y(:)) > R)
    error('Data out of range')
end

if any(sqrt(variance(:)) > sigma_max)
    error('sigma out of range: sigma_max = %f, max(variance) = %f\n', sigma_max, sqrt(max(variance(:))));
end

if length(variance) > n_s - 1
    error('too many components')
end

param_range = max(Y, [], 1) - min(Y, [], 1);
sigma = sqrt(variance);

val = quad(@coding_length, Delta_min, Delta_max, tol, 0, Y, n_s, sigma_max, nModes, param_range, sigma);


function val = coding_length(Deltas, Y, n_s, sigma_max, nModes, range, sigma)
%

% Do *NOT* change anything in here without profiling it before and
% afterwards.  Even apparantly innocuous changes such as re-arranging
% expressions or lines can have significant effects on this since it
% is called millions of times.  The JIT is very unpredictable and you
% will be wrong if you try to guess what it will do.  Insignificant
% changes here that only affect round-off error will change the
% results of the optimisation by suprisingly large amounts.

% The operands of all square roots and logarithms are positive
% numbers, but taking the real part reassures the JIT that the results
% are not complex.  Again: do *NOT* add or remove calls to `real'
% without profiling the code before and afterwards.

% This has been optimised for Matlab 6.5 on Linux/Intel.

% Distinguish between delta and Delta!  The best way is to make your editor
% display them as the appropriate Greek symbols.
nDeltas = length(Deltas);
val = zeros(nDeltas, 1);

ns_sqrt = real(sqrt(12/n_s));

for i=1:nDeltas
    Delta = Deltas(i);
    sigma_min = 2 * Delta;
    this_val = 0;
    
    for m=1:nModes,
        if range(m) < Delta
            % Case three
            continue
        end
        
        sigma_m = sigma(m);
        
        if sigma_m > sigma_min
            % Case one
            delta = min(1, sigma_m * ns_sqrt);
            
            L_data = -n_s * real(log(Delta)) ...
                + (n_s/2) * real(log(2*pi*sigma_m^2)) ...
                + (n_s/2) ...
                + (n_s * delta^2) / (12 * sigma_m^2);
        else
            % Case two
            
            % Quantize Ym.  Probably not necessary but we'll do it to be
            % on the safe side.
            Ym = Delta * round(Y(:, m) / Delta);
            delta = min(1, sigma_min * ns_sqrt);
            
            L_data = -n_s * real(log(Delta)) ...
                + (n_s/2) * real(log(2*pi*sigma_min^2)) ...
                + (1/(2 * sigma_min^2)) * sum(Ym.^2);
        end
        
        % $$$       L_delta = 1 + abs(real(log(delta)));
        % $$$       L_sigma = real(log((sigma_max - sigma_min) / delta));
        % $$$       L_parameters = L_sigma + L_delta;
        % $$$       L_total = L_parameters + L_data;
        this_val = this_val ...
            + real(log((sigma_max - sigma_min) / delta)) ...
            + 1 + abs(real(log(delta))) ...
            + L_data;
    end
    
    val(i) = this_val;
end
