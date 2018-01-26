function [] = write_tree_txt(tree, fname, out_dims)
%WRITE_CTREE_TXT write out a 2D array in ascii text format
%   [] = write_array_txt(complex_mat, fname)
%
% Inputs:
%      in_array - 2d array to be written
%
%      fname - filename to be written to
%
%
% Outputs: none
%
% Example:
%
% Notes: Each element 4 digits of
% precision. If, an element is less than 1e-4 it is written as 0 to save
% space
%
% See also:
%
% Created: 23-Apr-2015
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
if ~exist('out_dims', 'var')
    out_dims = [];
end

if ~isempty(out_dims)
    branches = tree.var > 0;
    
    if length(out_dims) == 1
        tree.var(branches) = tree.var(branches) + out_dims;
    else
         tree.var(branches) = out_dims(tree.var(branches));
    end
    tree.var(~branches) = -1;
end

cxx_names = {
    'node_features'
    'node_thresholds'
    'node_means'
    'node_vars'
    'node_children_left'
    'node_children_right'};
    
tree.children_l = tree.children(:,1);
tree.children_r = tree.children(:,2);

if strcmpi(tree.method, 'regression');
    fieldnames = {
        'var'
        'cut'
        'class'
        'nodeerr'
        'children_l'
        'children_r'};
    
    %Deal with complex valued orientation trees
    if any(~isreal(tree.class))
        tree.class = [real(tree.class) imag(tree.class)];
        tree.nodeerr = [tree.nodeerr zeros(size(tree.nodeerr))];
    end
    n_outputs = size(tree.class,2);
    
else %tree.method = 'classification'
    if size(tree.classprob, 2) ~= 2
        warning('Method only suitable for binary classsifcation trees');
    end
    
    fieldnames = {
        'var'
        'cut'
        'classprob'
        'nodeerr'
        'children_l'
        'children_r'};
    tree.classprob(:,end) = [];
    n_outputs = size(tree.classprob,2);    
end

formats = {'%d ', '%.4f ', '%.4f ', '%.4f ', '%d ', '%d'};  

n_features = tree.npred;
n_nodes = length(tree.node);
fid1 = fopen(fname, 'wt');

fprintf(fid1,'n_features: %d\n', n_features);
fprintf(fid1,'n_nodes: %d\n', n_nodes);
fprintf(fid1,'n_outputs: %d\n', n_outputs);

for i_f = 1:length(fieldnames)
    fprintf(fid1, '%s: {\n', cxx_names{i_f});

    for i_n = 1:n_nodes

        n_cols = size( tree.(fieldnames{i_f}), 2);
        
        for i_c = 1:n_cols
            fprintf(fid1, formats{i_f}, tree.(fieldnames{i_f})(i_n, i_c));
        end
        fprintf(fid1,'\n');
    end
    fprintf(fid1,'} \n');
end
fclose(fid1);




