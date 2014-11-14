function [forest_out] = mb_rf_disk_to_memory(forest_in)
%MB_RF_DISK_TO_MEMORY *Insert a one line summary here*
%   [forest_out] = mb_rf_disk_to_memory(forest_in)
%
% Inputs:
%      forest_in- *Insert description of input variable here*
%
%
% Outputs:
%      forest_out- *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 04-Feb-2010
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester

%Copy input forest structure to output forest structure and throw away
%unused fields
forest_out = rmfield(forest_in, {'tree_dir', 'tree_root'});

%For each tree
for ii = 1:length(forest_in.trees)
    
    %load tree and copy into output forest
    forest_out.trees{ii} = ...
        u_load([forest_in.tree_root forest_in.tree_dir forest_in.trees{ii}]);
    
end