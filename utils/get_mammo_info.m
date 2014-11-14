function [out1, out2, out3] = get_mammo_info(names, num_length)
%GET_MAMMO_INFO *Insert a one line summary here*
%   [out1, out2, out3] = get_mammo_info(num_length)
%
% Inputs:
%      names - *Insert description of input variable here*
%      
%      num_length - *Insert description of input variable here*
%
%
% Outputs:
%      out1 - *Insert description of input variable here*
%
%      out2 - *Insert description of input variable here*
%
%      out3 - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 22-Sep-2010
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester
if nargin < 2
    num_length = 3;
end

was_char = false;
if isstruct(names) && isfield(names, 'name')
    num_names = length(names);
    names_cell = cell(num_names,1);
    for ii = 1:num_names
        names_cell{ii} = names(ii).name;
    end
    names = names_cell; clear names_cell;
elseif ischar(names)
    num_names = 1;
    names = cellstr(names);
    was_char = true;
elseif iscellstr(names)
    num_names = length(names);
else
    error(['Data type of input names is not supported. Are you sure they are mammogram names?']);
end

nums = cell(num_names,1);
breasts = cell(num_names,1);
views = cell(num_names,1);

for ii = 1:num_names
    view_idx = strfind(names{ii}, 'CC');
    if isempty(view_idx)
        view_idx = strfind(names{ii}, 'ML');
    end
    if isempty(view_idx)
        error(['ML or CC not found in name ' num2str(ii) ': ' names{ii} '. Are you sure this is mammogram name?']);
    end
    nums{ii} = names{ii}(view_idx - 1 - num_length : view_idx - 2);
    breasts{ii} = names{ii}(view_idx - 1);
    views{ii} = names{ii}(view_idx:view_idx+1);
end
if was_char
    nums = char(nums);
    breasts = char(breasts);
    views = char(views);
end
if nargout == 1
    out1 = strcat(nums, breasts, views);
elseif nargout == 2
    out1 = nums;
    out2 = strcat(breasts, views);
elseif nargout == 3
    out1 = nums;
    out2 = breasts;
    out3 = views;
end