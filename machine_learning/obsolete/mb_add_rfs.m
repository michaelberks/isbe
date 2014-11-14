function [random_forest] = mb_add_rfs(source_rf_path, target_rf_path, delete_trees)
%MB_ADD_RFS *Insert a one line summary here*
%   [] = mb_combine_rfs(varargin)
%
% Inputs:
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 03-Feb-2010
% Author: Michael Berks
% Email : michael.berks@postgrad.man.ac.uk
% Phone : +44 (0)161 275 1241
% Copyright: (C) University of Manchester

% Set default for delete_trees
if nargin < 3
    delete_trees = 3;
end

%load in the source and target random forests
source_rf = u_load(source_rf_path);
target_rf = u_load(target_rf_path);

%otherwise check number of dimensions and classes are the same
if (source_rf.D ~= target_rf.D) || ...
        (isfield(target_rf, 'nclasses') && (source_rf.nclasses ~= target_rf.nclasses))
    %Warn user
    warning('Random forests do not appear to be compatible'); %#ok
    return;
end

curr_tree = length(target_rf.trees)+1;
    
%loop through each tree in loaded forest
for jj = 1:length(source_rf.trees)

    %save new tree name in main forest
    tree_name = ['rf_tree', zerostr(curr_tree, 4), '.mat'];
    target_rf.trees{curr_tree} = tree_name;

    %increment tree counter
    curr_tree = curr_tree + 1;

    %if asked, copy tree into new rf_directory
    copyfile(...
        [source_rf.tree_root source_rf.tree_dir source_rf.trees{jj}], ...
        [target_rf.tree_root target_rf.tree_dir tree_name]);

    %if asked delete the old tree
    if delete_trees
        delete([source_rf.tree_root source_rf.tree_dir source_rf.trees{jj}]);
    end
end
if delete_trees
    delete(source_rf_path);
    rmdir([source_rf.tree_root source_rf.tree_dir]);
end
random_forest = target_rf;
save(target_rf_path, 'random_forest');




