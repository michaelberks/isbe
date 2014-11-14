function [bin_counts] = mb_dual_tree_weighted_hist(data, bin_edges, phase_dims, mag_dims)
%MB_DUAL_TREE_WEIGHTED_HIST *Insert a one line summary here*
%   [bin_counts] = mb_dual_tree_weighted_hist(data,bin_centres)
%
% Inputs:
%      data- nx2L set of magnitude-(ILP)phase pair feature vectors to be
%            histogramed(?), where n is the number of feature points and L 
%            is the number of levels DT decomposition 
%           (i.e. we have phase and magnitude at each level, so dimension is 2L)  
%
%      bin_edges- 2 element cell array in which each element is either:
%                   i) a scalar specifying the number of bins between -pi/2,pi/2 or
%                   ii) a 1D vector specifying the bin edges
%                   along each phase dimension
%
%      phase_dims- 2 element vector specifying the two phase dimensions
%
%      mag_dims- 1D vector specifying the magnitudes to sum for a given
%                   phase-phase bin
%
%
% Outputs:
%      bin_counts- weighted-histogram counts for each phase-phase bin
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 20-Jan-2009
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester 

% Assume data is set of magnitude/phase pairs at each level of some
% dual-tree decomposition and thus mag/phase_dims should be doubled (e.g. level
% X phase is the in the 2*X column of data)
mag_dims = 2*mag_dims - 1;
phase_dims = 2*phase_dims; 

if length(bin_edges{1}) == 1
    num_bins1 = bin_edges{1};
    bin_edges{1} = linspace(-pi/2, pi/2, num_bins1+1);
else
    num_bins1 = length(bin_edges{1})-1;
end

if length(bin_edges{2}) == 1
    num_bins2 = bin_edges{2};
    bin_edges{2} = linspace(-pi/2, pi/2, num_bins2+1);
else
    num_bins2 = length(bin_edges{2})-1;
end

%Pre-allocate space for histogram counts
bin_counts = zeros(num_bins1, num_bins2);

% Compute histogram sums for all combinations of phase-phase
for p1 = 1:num_bins1
    for p2 = 1:num_bins2
        %Get indices for feature points included in this bin
        idx = ...
            data(:, phase_dims(1)) >= bin_edges{1}(p1) & ...
            data(:, phase_dims(1)) < bin_edges{1}(p1+1) & ...
            data(:, phase_dims(2)) >= bin_edges{2}(p2) & ...
            data(:, phase_dims(2)) < bin_edges{2}(p2+1);
        
        %Take sum of the nominated magnitudes for these points
        bin_counts(p1, p2) = sum(sum(data(idx, mag_dims)));
        
        %Below performs an unweighted count - this is obviously MUCH slower
        %than using Matlab's hist3, but is useful for checking we're
        %counting in the right bins
        %bin_counts(p1, p2) = sum(idx);
        
    end
end
        
        


