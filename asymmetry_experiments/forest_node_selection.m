function [dt_dims_counts dt_dims_scores dt_dims_cooc] = forest_node_selection(random_forest, make_figs, tree_root)
%FOREST_NODE_SELECTION *Insert a one line summary here*
%   [] = forest_node_selection(forest, make_figs)
%
% Inputs:
%      random_forest - *Insert description of input variable here*
%
%      make_figs - *Insert description of input variable here*
%
%
% Outputs:
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
if nargin < 2
    make_figs = 1;
end
if nargin < 3
    tree_root = [];
end
if ~isempty(tree_root)
    random_forest.tree_root = tree_root;
end

win_size = 9;
bands = 6;
num_levels = 4;
dt_dims_counts = zeros(2*bands*num_levels*win_size,1);
dt_dims_scores = zeros(2*bands*num_levels*win_size,1);
dt_dims_cooc = zeros(2*bands*num_levels*win_size);

%Loop through trees in forest
for jj = 1:length(random_forest.trees)
    
    %load tree
    tree = u_load([random_forest.tree_root random_forest.tree_dir random_forest.trees{jj}]);
    leaves = find(~any(tree.children,2));
    
    for ii = 1:length(leaves)
        %Get leaf node and score at this leaf
        curr_node = leaves(ii);
        score = abs(tree.class(curr_node));
        
        %For each leaf, follow the branch up through the tree to the root
        branch_dims = [];
        while(curr_node ~= 1)
            %Get the parent of the leaf and the splitting dimension of that
            %node
            curr_node = tree.parent(curr_node);
            %split_dim = ceil(tree.var(curr_node) / win_size);
            split_dim = tree.var(curr_node);
            
            %Increment the histogram counts and scores for this splitting
            %dimension
            dt_dims_cooc(branch_dims, split_dim) = ...
                dt_dims_cooc(branch_dims, split_dim) + 1;
            dt_dims_counts(split_dim) = dt_dims_counts(split_dim) + 1;
            dt_dims_scores(split_dim) = dt_dims_scores(split_dim) + score;
            
            %Add the splitting dimension to list of branch dimensions
            branch_dims(end+1) = split_dim; %#ok
        end
    end
end
%
if make_figs
    cmap = jet(9);

    f1 = figure; hold on; a1 = gca;
    f2 = figure; hold on; a2 = gca;
    f3 = figure; hold on; a3 = gca;
    for pos = 1:9
        idx = (0:9:431) + pos;
        bar(a1, idx, dt_dims_counts(idx), 'facecolor', cmap(pos,:), 'barwidth', 0.1);
        bar(a2, idx, dt_dims_scores(idx), 'facecolor', cmap(pos,:), 'barwidth', 0.1);
        bar(a3, idx, dt_dims_scores(idx)./dt_dims_counts(idx), 'facecolor', cmap(pos,:), 'barwidth', 0.1);
    end
    %
    y1 = get(a1, 'ylim');
    y2 = get(a2, 'ylim');
    y3 = get(a3, 'ylim');

    plot(a1, [(54.5:54:432)' (54.5:54:432)'], 0.95*y1, 'k:')
    plot(a2, [(54.5:54:432)' (54.5:54:432)'], 0.95*y2, 'k:')
    plot(a3, [(54.5:54:432)' (54.5:54:432)'], 0.95*y3, 'k:')

    plot(a1, [216.5 216.5], y1, 'r:')
    plot(a2, [216.5 216.5], y2, 'r:')
    plot(a3, [216.5 216.5], y3, 'r:')

    set(a1, 'xtick', 27:54.5:432, 'xticklabel', {'Level 1', 'Level 2', 'Level 3', 'Level 4'});
    axes(a1); text(100, y1(2), '\leftarrow Magnitude \rightarrow'); text(316, y1(2), '\leftarrow Phase \rightarrow');
    title(a1, 'Number of times each dimension is used as a splitting variable');
    set(a2, 'xtick', 27:54.5:432, 'xticklabel', {'Level 1', 'Level 2', 'Level 3', 'Level 4'});
    axes(a2); text(100, y2(2), '\leftarrow Magnitude \rightarrow'); text(316, y2(2), '\leftarrow Phase \rightarrow');
    title(a2, 'Sum of leaf dispersions for each dimension');
    set(a3, 'xtick', 27:54.5:432, 'xticklabel', {'Level 1', 'Level 2', 'Level 3', 'Level 4'});
    axes(a3); text(100, y3(2), '\leftarrow Magnitude \rightarrow'); text(316, y3(2), '\leftarrow Phase \rightarrow');
    title(a3, 'Average leaf dispersion for each dimension');
    
%     cmap = jet(6);
% 
%     f1 = figure; hold on; a1 = gca; %#ok
%     f2 = figure; hold on; a2 = gca; %#ok
%     f3 = figure; hold on; a3 = gca; %#ok
%     for band = 1:6
%         idx = (0:6:47) + band;
%         bar(a1, idx, dt_dims_counts(idx), 'facecolor', cmap(band,:), 'barwidth', 0.1);
%         bar(a2, idx, dt_dims_scores(idx), 'facecolor', cmap(band,:), 'barwidth', 0.1);
%         bar(a3, idx, dt_dims_scores(idx)./dt_dims_counts(idx), 'facecolor', cmap(band,:), 'barwidth', 0.1);
%     end
%     %
%     y1 = get(a1, 'ylim');
%     y2 = get(a2, 'ylim');
%     y3 = get(a3, 'ylim');
% 
% %     plot(a1, [(54.5:54:432)' (54.5:54:432)'], 0.95*y1, 'k:')
% %     plot(a2, [(54.5:54:432)' (54.5:54:432)'], 0.95*y2, 'k:')
% %     plot(a3, [(54.5:54:432)' (54.5:54:432)'], 0.95*y3, 'k:')
% % 
% %     plot(a1, [216.5 216.5], y1, 'r:')
% %     plot(a2, [216.5 216.5], y2, 'r:')
% %     plot(a3, [216.5 216.5], y3, 'r:')
% 
%     set(a1, 'xtick', 3.5:6:59, 'xticklabel', {'Level 1', 'Level 2', 'Level 3', 'Level 4'});
%     axes(a1); text(10, y1(2), '\leftarrow Magnitude \rightarrow'); text(40, y1(2), '\leftarrow Phase \rightarrow');
%     title(a1, 'Number of times each dimension is used as a splitting variable');
%     set(a2, 'xtick', 3.5:6:59, 'xticklabel', {'Level 1', 'Level 2', 'Level 3', 'Level 4'});
%     axes(a2); text(10, y2(2), '\leftarrow Magnitude \rightarrow'); text(40, y2(2), '\leftarrow Phase \rightarrow');
%     title(a2, 'Sum of leaf dispersions for each dimension');
%     set(a3, 'xtick', 3.5:6:59, 'xticklabel', {'Level 1', 'Level 2', 'Level 3', 'Level 4'});
%     axes(a3); text(10, y3(2), '\leftarrow Magnitude \rightarrow'); text(40, y3(2), '\leftarrow Phase \rightarrow');
%     title(a3, 'Average leaf dispersion for each dimension');
end
