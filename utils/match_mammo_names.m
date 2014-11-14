function [matched_names, missing_idx] = match_mammo_names(in_dir, mammo_names, post_str, pre_str, file_ext)
%MATCH_MAMMO_NAMES *Insert a one line summary here*
%   [matched_names, idx] = match_mammo_names(in_dir, mammo_names, post_str, pre_str, file_ext)
%
% Inputs:
%      in_dir - *Insert description of input variable here*
%
%      mammo_names - *Insert description of input variable here*
%
%      post_str - *Insert description of input variable here*
%
%      pre_str - *Insert description of input variable here*
%
%      file_ext - *Insert description of input variable here*
%
%
% Outputs:
%      matched_names - *Insert description of input variable here*
%
%      idx - *Insert description of input variable here*
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

if nargin < 5
    file_ext = '.mat';
end
if nargin < 4
    pre_str = [];
end
if nargin < 3
    post_str = [];
end
if nargin < 2
    error('At least 2 arguments (in_dir and mammo_names) must be specified');
end

%Get full list of files from input directory
full_list = dir([in_dir '/*' pre_str '*' post_str file_ext]);

%Get mammo names from full list
full_names = get_mammo_info(full_list);

%Get the indices of full names that match the input mammo names
[tf loc] = ismember(full_names, mammo_names);

%Create output cell of matched names
matched_names = cell(length(mammo_names), 1);
for ii = 1:length(tf)
    if tf(ii)
        matched_names{loc(ii)} = full_list(ii).name;
    end
end

%Get missing idx
missing_idx = setdiff(1:length(mammo_names), loc(tf))';
