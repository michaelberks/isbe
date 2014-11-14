function [full_ilp] = ilp_same_size(full_tree)
%ILP_SAME_SIZE take inter-level product of dual-tree that has already been
%interpolated to the full pixel grid grid
%   [full_ilp] = ilp_same_size(full_tree)
%
% Inputs:
%      full_tree- dual-tree that has already been
%       interpolated to the full pixel grid grid - (m,n,6,l) array
%
%
% Outputs:
%      full_ilp- (m,n,6,L) array, such that the top L-1 levels are products
%        of the original full_tree level l, and the conjugate of phase-doubled
%        level l+1
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 11-May-2009
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester 

full_ilp = full_tree; clear full_tree;


%for each level from fine to coarse subtract the
%phase doubled coarser level - note we divide by magnitudes of the coarser
%level so the resulting product has the same magnitude as the finer level
for l = 1:size(full_ilp,4)-1
    
    full_ilp(:,:,:,l) = full_ilp(:,:,:,l) .* conj(full_ilp(:,:,:,l+1).^2) ...
                            ./ (abs(full_ilp(:,:,:,l+1)).^2 + 1e-6);
end