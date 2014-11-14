function [matched_names1, matched_names2 matched_names] = match_generic_names(file_list1, file_list2)
%MATCH_MAMMO_NAMES *Insert a one line summary here*
%   [matched_names, idx] = match_file_names(in_dir, file_names, post_str, pre_str, file_ext)
%
% Inputs:
%      file_list1 - *Insert description of input variable here*
%
%      file_list2 - *Insert description of input variable here*
%
%
% Outputs:
%      matched_names1 - *Insert description of input variable here*
%
%      matched_names2 - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 16-Nov-2010
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester

%Get names from two lists
[file_names1 pre1 post1] = get_generic_names(file_list1);
[file_names2 pre2 post2] = get_generic_names(file_list2);

%Get the indices of full names that match the input mammo names
[matched_names] = intersect(file_names1, file_names2);

%Create output cells of matched names
matched_names1 = strcat(pre1, matched_names, post1);
matched_names2 = strcat(pre2, matched_names, post2);
