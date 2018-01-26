function [] = write_forest_txt(forest, out_dir, out_dims)
%WRITE_FOREST_TXT write out a decision tree in text format
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

out_dir = [out_dir '\'];
create_folder(out_dir);

fname = [out_dir 'tree_list.txt'];
fid = fopen(fname, 'wt');

for i_tree = 1:length(forest.trees)
    if isstruct(forest.trees{i_tree})
        tree = forest.trees{i_tree};
        tree_name_out = [out_dir 'rf_tree' zerostr(i_tree,4) '.txt'];      
    else
        tree_name_in = [forest.tree_root forest.tree_dir forest.trees{i_tree}];
        tree = u_load(tree_name_in);
        tree_name_out = [out_dir forest.trees{i_tree}(1:end-3) 'txt'];
    end
      
    write_tree_txt(tree, tree_name_out, out_dims);
    
    find_backslash = tree_name_out=='\';
    tree_name_out(find_backslash) = '/';
    
    fprintf(fid, '%s', tree_name_out);
    fprintf(fid,'\n');
end
fclose(fid);