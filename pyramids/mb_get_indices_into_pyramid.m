function [r, c] = mb_get_indices_into_pyramid(varargin)
%
% [r, c] = cjr_get_indices_into_pyramid(...)
%
% Computes indices into pyramids.
%
% This function operates in two modes to achieve two different
%	aims.
%
% Mode 1:
% -------
%
% [r, c] = cjr_get_indices_into_pyramid('Function', 'UpPyramid', ...
%						'Row', r, 'Col', c, 'Level', l)
%
% Finds the indices into a particular level that correspond to indices
% into the finest level of a pyramid (a quadtree).
%
% r and c are row and column indices which start at 1.
%
% l is the pyramid level that you want the corresponding
% indices for. The bottom level of the pyramid is level 1 and
% is the finest level.
%
% r and c may be column vectors if you want to calculate multiple
% indices.
%
% Mode 2:
% -------
%
% [r, c] = cjr_get_indices_into_pyramid('Function', 'DownPyramid', ...
%						'Row', r, 'Col', c, 'LevelsBelow', l)
%
% Finds indices l levels below the level containing the indices (r, c)
% that are on a path through (r, c) to the top of the pyramid. Note
% that for each element in a pyramid level, there are four child nodes,
% and so it is IMPOSSIBLE to uniquely determine paths from higher levels to lower
% levels by working down the pyramid!
%
% CJR_SAMPLE_SYN_MAMMO_TEXTURE uses the U_PACKARGS
% interface function and so arguments to the function may be passed
% as name-value pairs or as a struct with fields with corresponding
% names. The field names are defined as:
%
% Mandatory Arguments:
%
%	'Function'
%		- determines whether to find the indices into a particular level that
%		correspond to indices into the finest level of a pyramid ('UpPyramid'),
%		or whether to find indices at the base of the pyramid that are on a path
%		from a pair of given indices at a particular level ('DownPyramid'). Specify either:
%			'UpPyramid', or
%			'DownPyramid'
%
% 'Row'
%		- the row (what this memans depends upon the mode used).
%
% 'Col'
%		- the column (what this memans depends upon the mode used).
%
% Other arguments (see documentation above to determine which ones to use):
%
% 'Level'
%		- in mode 1, this is the pyramid level that you want the corresponding
% 	indices for. The bottom level of the pyramid is level 1 and
% 	is the finest level.
%
%	'LevelsBelow'
%		- used by mode 2. See the description above. An ASCII diagram might help:
%
%                   ()
% ->              ()  (X)
%             ()  ()  (Y)  (Y)
%     ()  ()  ()  ()  (Y)  (Y)  (Y)  (Y)
%
% The above diagrams is a bi-tree pyramid (i.e. 1D images at each level). Imagine we have an
% index for the level indicated by the -> arrow; the index has the value 2 (indicated y the
% position of the X). We want to know an index value for the bottom level of the pyramid that
% is on the path through our index to the top of the pyramid. There are several such paths,
% but we only want an index into the bottom level that specifies one such path. Mode 2 of this
% function returns such an index. The Ys indicate nodes that solve the above problem, and
% demonstrate that it is IMPOSSIBLE to uniquely determine paths from higher levels to lower
% levels by working down the pyramid.
%
% See also: CJR_PYR_PARENT_INDICES

% unpack the arguments
args = u_packargs(varargin, 'strict', ...
			{'Function', ... % Mandatory arguments
			'Row', ...
			'Col'}, ... 
			'Level',  [], ... % Optional arguments
			'LevelsBelow', [] ...
			);

if strcmp(upper(args.Function), 'UPPYRAMID')
	if args.Level == 1
        r = args.Row;
        c = args.Col;
    else
        r = ceil(args.Row/(2^(args.Level-2)));
        c = ceil(args.Col/(2^(args.Level-2)));
    end
elseif strcmp(upper(args.Function), 'DOWNPYRAMID')
	r = 2 * args.LevelsBelow * args.Row;
	c = 2 * args.LevelsBelow * args.Col;
else
	error(['Mode named ' args.Function ' is not supported']);
end