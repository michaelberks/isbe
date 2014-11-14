function [leaf_output] = forest_leaf_output(random_forest, varargin)
%FOREST_LEAF_OUTPUT *Insert a one line summary here*
%   [leaf_hist] = forest_leaf_output(forest, orientation, make_figs)
%
% Inputs:
%      forest - *Insert description of input variable here*
%
%      make_figs - *Insert description of input variable here*
%
%
% Outputs:
%      leaf_hist - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 15-Feb-2011
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
args = u_packargs(varargin,... % the user's input
         '0', ... % non-strict mode
         'tree_root', [],...
         'make_figs', 0,... % the optional arguments
         'do_dang', 1,...
         'd_bins', linspace(0,1,50),...
         'size_bins', linspace(2,200,100));
clear varargin;        

if ~isempty(args.tree_root)
    random_forest.tree_root = args.tree_root;
end

num_trees = length(random_forest.trees);

%pre-allocate histogram counts
d_hist = zeros(length(args.d_bins),1);
size_hist = zeros(length(args.size_bins),1);

ang_hist = zeros(180,1);
dang_hist = zeros(180,1);

tree_mean_d = zeros(num_trees,1);
tree_median_d = zeros(num_trees,1);
tree_levels = zeros(num_trees,1);
tree_sizes = zeros(num_trees,1);
%Loop through trees in forest
for ii = 1:num_trees
    
    %load tree
    tree = u_load([random_forest.tree_root random_forest.tree_dir random_forest.trees{ii}]);
    
    %Get number of levels in tree
    tree_sizes(ii) = length(tree.node);
    tree_level = nan(tree_sizes(ii),1);
    tree_level(1) = 0;
    for jj = 2:tree_sizes(ii)
        tree_level(jj) = tree_level(tree.parent(jj))+1;
    end
    tree_levels(ii) = max(tree_level);
    
    %find leaves
    leaves = ~tree.var;
    scores = abs(tree.class(leaves));
    
    tree_mean_d(ii) = mean(abs(tree.class(leaves)));
    tree_median_d(ii) = median(abs(tree.class(leaves)));
    
    %find the output of each leaf to the nearest degree and increment the
    %appropriate histogram count
    if args.do_dang
        degs = mod(round(90*angle(tree.class(leaves))/pi)-1,180)+1;

        for jj = 1:length(degs)
            ang_hist(degs(jj)) = ang_hist(degs(jj)) + 1;
            dang_hist(degs(jj)) = dang_hist(degs(jj)) + scores(jj);
        end
    else
        ang_hist_tree = hist(mod(90*angle(tree.class(leaves))/pi,180), 1:180);
        ang_hist = ang_hist + ang_hist_tree(:);
    end    
    
    %Increment histogram counts for dispersion and leaf size
    d_hist_tree = hist(scores, args.d_bins);
    d_hist = d_hist + d_hist_tree(:);
    
   size_hist_tree = hist(tree.nodesize(leaves), args.size_bins);
   size_hist = size_hist + size_hist_tree(:);
end

leaf_output.ang_hist = ang_hist;
leaf_output.d_hist = d_hist;
leaf_output.dang_hist = dang_hist;
leaf_output.size_hist = size_hist;
leaf_output.tree_mean_d = tree_mean_d;
leaf_output.tree_median_d = tree_median_d;
leaf_output.tree_levels = tree_levels;
leaf_output.tree_sizes = tree_sizes;

if args.make_figs
%     figure; bar(1:180, leaf_hist);
%     title('Histogram of forest leaf outputs')
%     xlabel('Leaf output (degrees)');
end

