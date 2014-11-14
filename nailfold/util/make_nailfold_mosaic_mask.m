function [nailfold_mask] = make_nailfold_mosaic_mask(nailfold, g_lim, border_sz)
%MAKE_NAILFOLD_MOSAIC_MASK *Insert a one line summary here*
%   [nailfold_mask] = make_nailfold_mosaic_mask(nailfold, g_lims)
%
% Inputs:
%      nailfold - *Insert description of input variable here*
%
%      g_lim - *Insert description of input variable here*
%
%
% Outputs:
%      nailfold_mask - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 12-Mar-2013
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
if nargin < 2
    g_lim = 250;
end
if nargin < 3
    border_sz = 10;
end

nailfold_mask = nailfold > g_lim;
edges = false(size(nailfold));
edges([1 end], :) = 1;
edges(:, [1 end]) = 1;
[edge_r edge_c] = find(edges);

nailfold_mask = ~bwselect(nailfold_mask, edge_c, edge_r);

if border_sz
    nailfold_mask = imerode(nailfold_mask, strel('disk', border_sz));
end
