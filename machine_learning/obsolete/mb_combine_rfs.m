function [random_forest] = mb_combine_rfs(varargin)
%MB_COMBINE_RFS *Insert a one line summary here*
%   [] = mb_combine_rfs(varargin)
%
% MB_COMBINE_RFS uses the U_PACKARGS interface function
% and so arguments to the function may be passed as name-value pairs
% or as a struct with fields with corresponding names. The field names
% are defined as:
%
% Mandatory Arguments:
%
% Optional Arguments:
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

% Unpack the arguments:
args = u_packargs(varargin,... % the user's input
    '0', ... % non-strict mode
    {'rf_dir'},...% the mandatory arguments
    'tree_dir', [],...
    'tree_root',[],...
    'rf_list', [],...
    'replace_tree_root', [],...
    'copy_trees', 1,...% the optional arguments
    'delete_trees', 0,...
    'save_path', []);
clear varargin;

%if we're going to copy trees we need to make sure we have directory to put
%them in
if args.copy_trees
    tree_dir = [args.tree_root args.tree_dir];
    
    if isempty(tree_dir)
        warning(['You have selected to copy trees but not specified an output directory.',...
            ' New trees will be copied to the current directory']); %#ok
    else
        %make sure tree dir is filesep terminated
        if ~strcmp(tree_dir(end), filesep)
            tree_dir = [tree_dir filesep];
        end
        %if it isn't already a directory create it
        if ~exist(tree_dir, 'dir')
            mkdir(tree_dir);
			if ~ispc
				fileattrib(tree_dir,'+w','g');
			end
        end
    end
    random_forest.tree_dir = args.tree_dir;
    random_forest.tree_root = args.tree_root;
    
elseif args.delete_trees
    %if we're not copying trees, make sure we don't delete them
    warning('You should not delete trees you haven''t copied. Turning delete_trees off'); %#ok
    args.delete_trees = 0;
end

%make sure rf_dir is filesep terminated
rf_dir = args.rf_dir;
if isempty(rf_dir) && ~strcmp(rf_dir(end), filesep)
    rf_dir = [rf_dir filesep];
end

%If we haven't been given a list of forests to combine, get directory
%listing of forests
rf_list = args.rf_list;
if isempty(rf_list)
    rf_list = dir([rf_dir '*random_forest*.mat']);
end

n_rfs = length(rf_list);
curr_tree = 1;

%go through each forest adding each tree to main forest
for ii = 1:n_rfs
    
    rf = u_load([rf_dir rf_list(ii).name]);
    
    %if first forest get information and save in main forest
    if ii == 1
        random_forest.D = rf.D;
        random_forest.d = rf.d;
        
        %These fields are only for classification forests
        if isfield(rf, 'nclasses')
            random_forest.nclasses = rf.nclasses;
            random_forest.classname = rf.classname;
            random_forest.prior = rf.prior;
            random_forest.cost = rf.cost;
            random_forest.names = rf.names;
        end
        
        
        %otherwise check number of dimensions and classes are the same
    elseif (random_forest.D ~= rf.D) || ...
            (isfield(random_forest, 'nclasses') && (random_forest.nclasses ~= rf.nclasses))
        %Warn user
        warning('Random forests do not appear to be compatible'); %#ok
        display(['Skipping forest ', num2str(ii), ': ', rf_list(ii).name]);
        continue;
    end
    
    %check if we need to replace the root dircetory of the tree (a bit of a
    %fudge to make it easier to copy trees built on a network cluster
    %machine)
    if args.copy_trees && ~isempty(args.replace_tree_root)
        rf.tree_root = args.replace_tree_root;
    end
    
    %loop through each tree in loaded forest
    for jj = 1:length(rf.trees)
        
        %save new tree name in main forest
        tree_name = ['rf_tree', zerostr(curr_tree, 4), '.mat'];
        random_forest.trees{curr_tree} = tree_name;
        
        %increment tree counter
        curr_tree = curr_tree + 1;
        
        %if asked, copy tree into new rf_directory
        if args.copy_trees
            copyfile([rf.tree_root rf.tree_dir rf.trees{jj}],[tree_dir tree_name]);
        end
        
        %if asked delete the old tree
        if args.delete_trees
            delete([rf.tree_root rf.tree_dir rf.trees{jj}]);
        end
    end
    if args.delete_trees
        delete([rf_dir rf_list(ii).name]);
        rmdir([rf.tree_root rf.tree_dir], 's');
    end
end

save_path = args.save_path;
if isempty(save_path)
    warning(['You have selected to save random_forest but not specified an output directory.',...
        ' random_forest will be save to the current directory']); %#ok
    save('random_forest', 'random_forest');
else
    save([save_path, 'random_forest.mat'], 'random_forest');
end




