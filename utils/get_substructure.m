function [str_out] = get_substructure(str_in, fields, str_out)
%GET_SUBSTRUCTURE Copies fields listed in inptu variable 'fields' from
%'str_in' to the returned 'str_out'
%   [str_out] = get_substructure(str_in, fields)
%
% Inputs:
%      str_in - Input structure
%
%      fields - Fields you want to keep
%
%
% Outputs:
%      str_out - Output structure with retained fields
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 03-Feb-2011
% Author: Phil Tresadern 
% Email : philip.tresadern@manchester.ac.uk 
% Phone : +44 (0)161 275 5114 
% Copyright: (C) University of Manchester 

% If str_out undefined then start from an empty structure
if ~exist('str_out','var'), str_out = []; end

for f = 1:length(fields)
    fieldname = fields{f};
	if isfield(str_in, fieldname)
		str_out.(fieldname) = str_in.(fieldname);
	else
		warning('ASYM:UnknownField',...
				['Field not found: ', fieldname]);
	end
end
