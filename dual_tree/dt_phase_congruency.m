function [phase_map] = dt_phase_congruency(image_in, min_level, max_level)
%DT_PHASE_CONGRUENCY *Insert a one line summary here*
%   [phase_map] = dt_phase_congruency(image_in,levels)
%
% Inputs:
%      image_in- *Insert description of input variable here*
%
%      levels- *Insert description of input variable here*
%
%
% Outputs:
%      phase_map- *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 07-May-2009
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester

%build dual-tree
dual_tree = dtwavexfm2(image_in, max_level+1);

%convert to the full tree
full_tree = dt_to_full_image(dual_tree, min_level:max_level+1);
clear dual_tree;

%Do the ILP stuff - i.e. for each level from coarse to fine subtract the
%phase doubled coarser level
for l = max_level-min_level+1:-1:1
    
    full_tree(:,:,:,l) = full_tree(:,:,:,l) .* conj(full_tree(:,:,:,l+1).^2) ...
                            ./ (abs(full_tree(:,:,:,l+1)).^2 + 1e-6);
end

%Discard the coarsest level and take the sum
full_tree(:,:,:,end) = [];
phase_map = sum(full_tree, 4);
        
    
    
