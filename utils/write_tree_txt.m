function [] = write_tree_txt(tree, fname)
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
    
     
end

formats = {'%d ', '%.4f ', '%.4f ', '%.4f ', '%d ', '%d'};  

n_features = tree.npred;
n_nodes = length(tree.node);

fid1 = fopen(fname, 'wt');

fprintf(fid1,'n_features: %d\n', n_features);
fprintf(fid1,'n_nodes: %d\n', n_nodes);

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




