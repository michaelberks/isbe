function [depth] = get_tree_depth(tree)
%GET_TREE_DEPTH *Insert a one line summary here*
%   [depth] = get_tree_depth(tree)
%
% Inputs:
%      tree - *Insert description of input variable here*
%
%
% Outputs:
%      depth - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 31-Aug-2010
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester 

%Get the leaves of the tree
leaves = find(~tree.children(:,1));

%For each leave work out the number of levels to the parent
depth = 0;
for ii = 1:length(leaves);

    parent = tree.parent(leaves(ii));
    count = 1;
    
    while parent
        count = count+1;
        child = parent;
        parent = tree.parent(child);
    end
    
    %If number of levels for this leaf is maximum save as depth variable
    depth = max(depth, count);
end