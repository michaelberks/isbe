function [tree_args] = get_tree_args_from(args_in, tree_args)
% Given a set of arguments, return only those arguments that are 
% relevant to decision trees (or append them if tree_args is supplied)

% deal with undefined inputs
if ~exist('tree_args','var'), tree_args = []; end

tree_args.random_m = args_in.d;
tree_args.prior = args_in.prior;
tree_args.cost = args_in.cost;
tree_args.split_criterion = args_in.split_criterion;
tree_args.split_min = args_in.split_min;
tree_args.prune = args_in.prune;
tree_args.names = args_in.names;