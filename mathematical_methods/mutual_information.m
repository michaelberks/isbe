function [MI, hist_xy] = mutual_information(x, y, n_bins_x, n_bins_y)
%MUTUAL_INFORMATION *Insert a one line summary here*
%   [MI] = mutual_information(x, y, n_bins_x, n_bins_y)
%
% Inputs:
%      x - *Insert description of input variable here*
%
%      y - *Insert description of input variable here*
%
%      n_bins_x - *Insert description of input variable here*
%
%      n_bins_y - *Insert description of input variable here*
%
%
% Outputs:
%      MI - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 25-Apr-2018
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester
if exist('n_bins_x', 'var') && ~isempty(n_bins_x)
    if exist('n_bins_y', 'var') && ~isempty(n_bins_y)
        hist_xy = histcounts2(x, y, [n_bins_x n_bins_y]);
    else
        hist_xy = histcounts2(x, y, [n_bins_x n_bins_x]);
    end
else
    hist_xy = histcounts2(x, y);
end
pxy = hist_xy / sum(hist_xy(:));

px = sum(pxy, 2); % marginal for x over y
py = sum(pxy, 1); % marginal for y over x
px_py = px .* py; % Broadcast to multiply marginals
% Now we can do the calculation using the pxy, px_py 2D arrays
nzs = pxy > 0;% # Only non-zero pxy values contribute to the sum
MI = sum(pxy(nzs) .* log(pxy(nzs) ./ px_py(nzs)));
