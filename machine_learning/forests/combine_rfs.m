function [random_forest] = combine_rfs(varargin)
%COMBINE_RFS *Insert a one line summary here*
%   [] = combine_rfs(varargin)
%
% COMBINE_RFS uses the U_PACKARGS interface function
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
    'sampled_maps_dir', [], ...
    'rf_list', [],...
    'rf_name', 'predictor',...
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
        %make sure tree dir is properly formed
        tree_dir = prettypath(tree_dir);
        
        %if it isn't already a directory create it
        create_folder(tree_dir);
    end
    random_forest.tree_dir = args.tree_dir;
    random_forest.tree_root = args.tree_root;
    
elseif args.delete_trees
    %if we're not copying trees, make sure we don't delete them
    warning('You should not delete trees you haven''t copied. Turning delete_trees off'); %#ok
    args.delete_trees = 0;
end

%Create directory for the sampled maps
if ~isempty(args.sampled_maps_dir)
    sampled_maps_dir = [args.tree_root args.sampled_maps_dir];
    create_folder(sampled_maps_dir);
    sampled_pts_list = [];
else
    sampled_maps_dir = [];
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
    rf_list = dir([rf_dir '*' args.rf_name '*.mat']);
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
        random_forest.sampled_data_dir = sampled_maps_dir;
        
        %These fields are only for classification forests
        if isfield(rf, 'nclasses')
            random_forest.nclasses = rf.nclasses;
            random_forest.classname = rf.classname;
            random_forest.prior = rf.prior;
            random_forest.cost = rf.cost;
            random_forest.names = rf.names;
        else
            random_forest.regression_method = rf.regression_method;
        end
        
        %If we're making sample maps copy the image list
        if ~isempty(sampled_maps_dir)
            image_lists = u_load([rf.tree_root rf.sampled_data_dir 'image_lists.mat']);
            [image_list{1:length(image_lists),1}] = image_lists.image; %#ok
            clear image_lists;
        end
        
        %otherwise check number of dimensions and classes are the same
    elseif (random_forest.D ~= rf.D) || ...
            (isfield(random_forest, 'nclasses') && (random_forest.nclasses ~= rf.nclasses))
        %Warn user
        warning('Random forests do not appear to be compatible'); %#ok
        display(['Skipping forest ', num2str(ii), ': ', rf_list(ii).name]);
        continue;
    end
    
    %Compute number of trees in this forest
    n_trees = length(rf.trees);
    
    %Get a list of sampled data files for this forest
    if ~isempty(sampled_maps_dir)
        sampled_pts_files = dir([rf.sampled_data_dir 'sampled_pts*.mat']);
        if length(sampled_pts_files) ~= n_trees
            display(['Warning, incorrect number of sampled pts files found',...
                'combine function will proceed without making sample maps']);
            sampled_maps_dir = [];
        end
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
        tree_name = ['rf_tree' zerostr(curr_tree, 4), '.mat'];
        random_forest.trees{curr_tree} = tree_name;
        
        %if asked, copy tree into new rf_directory
        if args.copy_trees
            copyfile([rf.tree_root rf.tree_dir rf.trees{jj}],[tree_dir tree_name]);
        end
        
        %if asked delete the old tree
        if args.delete_trees
            delete([rf.tree_root rf.tree_dir rf.trees{jj}]);
        end
        
        %Move sampled data files
        if ~isempty(sampled_maps_dir)
            sampled_pts_name = ['sampled_pts' zerostr(curr_tree, 4), '.mat'];
            sampled_pts_list{curr_tree} = [sampled_maps_dir sampled_pts_name]; %#ok
            if args.delete_trees
                movefile(...
                    [rf.tree_root rf.sampled_data_dir sampled_pts_files(jj).name],...
                    sampled_pts_list{curr_tree});
            else
                copyfile(...
                    [rf.tree_root rf.sampled_data_dir sampled_pts_files(jj).name],...
                    sampled_pts_list{curr_tree});
            end
        end
            
        %increment tree counter
        curr_tree = curr_tree + 1;
                
    end
    if args.delete_trees
        %Delete the forest then remove the tree dir
        delete([rf_dir rf_list(ii).name]);
        rmdir([rf.tree_root rf.tree_dir], 's');
        
        if ~isempty(sampled_maps_dir)
            %Remove the imae lists file and delete the old sampled data dir
            %(which should now be empty)
            delete([rf.tree_root rf.sampled_data_dir 'image_lists.mat']);
            rmdir([rf.tree_root rf.sampled_data_dir], 's');
        end
    end
end

if ~isempty(sampled_maps_dir)
    %Finally we can make the sample maps from the sampled pts files
    make_sampled_maps(image_list, sampled_pts_list, sampled_maps_dir, 1);
end

save_path = args.save_path;
if isempty(save_path)
    warning(['You have selected to save random_forest but not specified an output directory.',...
        ' random_forest will be save to the current directory']); %#ok
    save('random_forest', 'random_forest');
else
    save(save_path, 'random_forest');
end




