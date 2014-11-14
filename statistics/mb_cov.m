function C = mb_cov(X)
%MB_COV compute covariance of a matrix
%   [C] = mb_cov(X)
%
% Inputs:
%      X - input matrix
%
%      C - Output covariance matrix
%
%
% Outputs:
%      probability_image - Pixel-wise probability of belonging to a bar
%
%
% Example:
%
% Notes: If the number of data points in X is small, returns the standard
% non-biased covariance matrix. However if N is large...
%
% See also:
%
% Created: 27-Jan-2010
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester

%get size of X    
[n d] = size(X);

%if number of samples in X is less than 1e5 return standard unbiased
%covariance (i.e. divide Matlab's cov by (n-1)
if n < 1e5
    C = cov(X)*(n-1) / n;
else
%However, if n large...
    C = zeros(d);
    X_mean = mean(X);

    for i = 1:d
        for j= 1:i

            C(i,j) = (X(:,i)'-X_mean(i)) * (X(:,j) - X_mean(j));
            C(j,i) = C(i,j);

        end
    end

    C = C / n;

end
        
        