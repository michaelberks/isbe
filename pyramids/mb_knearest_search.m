function [values, idx, knearest] = mb_knearest_search(varargin)
%
% MB_KNEAREST_SEARCH Synthesise 
%
% MB_KNEAREST_SEARCH uses the U_PACKARGS interface function
% and so arguments to the function may be passed as name-value pairs
% or as a struct with fields with corresponding names. The field names
% are defined as:
%
% Mandatory Arguments:
%
% Return Value:
%
%   MB_KNEAREST_SEARCH returns ??.
%
% Example Usage:
%
% References:
%
% See also: MB_KNEAREST_BUILD

% Unpack the arguments:
args = u_packargs(varargin, '0', ... % strict mode
		  {...  % The mandatory arguments
          'Data',...
          'TestVector'
          },... % The optional arguments
          'Method', 'standard'...
          );

clear varargin;

data_length = size(args.Data.Data, 1);

switch args.Method
    case 'standard'
        
        dist = sum((args.Data.Data - repmat(args.TestVector, data_length, 1)).^2, 2);
        [dummy, idx] = min(dist);
        knearest = args.Data.Data(idx, :);
        
    case 'tree'

        [knearest idx] = mb_tsvq_search('Tree', args.Data.Data,...
            'TreeIndices', args.Data.TreeIndices,...
            'TestVector', args.TestVector);
        
end
% comparison = [args.TestVector' knearest']
values = args.Data.Values(idx,:);




        
        
        
        
        
        
        
        
        
        
        
        
        
        
        