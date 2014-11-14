function [codeword, idx] = mb_tsvq_search(varargin)

%MB_TSVQ_SEARCH Search a pre-build TSVQ tree
%
%
% MB_TSVQ_SEARCH uses the U_PACKARGS interface function
% and so arguments to the function may be passed as name-value pairs
% or as a struct with fields with corresponding names. The field names
% are defined as:
%
% Mandatory Arguments:
%	'Tree' - the TSVQ to search
%
%   'Test vector' - the test vector we're classifying
%
% Optional Arguments:
%	
%	BackTrack - Allow tree backtracking - not yet functional
%   
% MB_K_MEANS_CLUSTERING returns:
%	'codeword'
%			- The closest match in the tree 
%
%   'idx' - The index (possibly indices) of the codeword in the original data set
% References:
%
% pack the args
args = u_packargs(varargin,...
				'0',...
				{'Tree',... % the mandatory arguments
                 'TreeIndices',...
                 'TestVector',...
                 }, ... 
				'BackTrack', false);
clear varargin;

d1 = sum((args.TestVector - args.Tree(1, 4:end)).^2);
d2 = sum((args.TestVector - args.Tree(2, 4:end)).^2);

if d1 < d2
    next_node = 1;
else
    next_node = 2;
end
    
while true;
    this_node = next_node;
    
    if args.Tree(this_node, 1) < 0;
        break; %we've reached a leaf
    end
    
    d1 = sum((args.TestVector - args.Tree(args.Tree(this_node,1), 4:end)).^2);
    d2 = sum((args.TestVector - args.Tree(args.Tree(this_node,2), 4:end)).^2);
     
    if d1 < d2
        next_node = args.Tree(this_node, 1);
    else
        next_node = args.Tree(this_node, 2);
    end
end
codeword = args.Tree(this_node, 4:end);

idx = args.TreeIndices(this_node).Indices;
if length(idx) > 1
    idx = randsample(idx, 1);
end