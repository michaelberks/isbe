function [names_out pre_const post_const] = get_generic_names(names)
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
if isstruct(names) && isfield(names, 'name')
    num_names = length(names);
    names_cell = cell(num_names,1);
    for ii = 1:num_names
        names_cell{ii} = names(ii).name;
    end
    names = char(names_cell);
     
elseif iscellstr(names)
    names = char(names);
elseif ~ischar(names)
    error('Data type of input names is not supported. Are you sure they are image names?');
end

ascii_names = double(names);
constant_columns = ~any(diff(ascii_names));

pre_const = [];
if nargout > 1 && constant_columns(1)
    pre_col = find(~constant_columns, 1)-1;
    pre_const = cellstr(names(1,1:pre_col));
end
post_const = [];
if nargout > 2 && constant_columns(end)
    post_col = find(~constant_columns, 1, 'last')+1;
    post_const = cellstr(names(1,post_col:end));
end

names(:,constant_columns) = [];
names_out = cellstr(names);

