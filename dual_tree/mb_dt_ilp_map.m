function [ilp_map] = mb_dt_ilp_map(varargin)
%MB_DT_ORIENTATION_MAP *Insert a one line summary here*
%   [orientation,orientation_levels] = mb_dt_orientation_sum(varargin)
%
% MB_DT_ORIENTATION_MAP uses the U_PACKARGS interface function
% and so arguments to the function may be passed as name-value pairs
% or as a struct with fields with corresponding names. The field names
% are defined as:
%
% Mandatory Arguments:
%
% Optional Arguments:
%
% Outputs:
%      orientation- *Insert description of input variable here*
%
%      orientation_levels- *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 27-Mar-2009
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester 

% Unpack the arguments:
args = u_packargs(varargin, '0', ...
    'DualTree', [],...
    'Image', [],...
    'Subs', [],...
    'Levels', 2:5);

%Check the user has supplied either an image or a dual_tree
if isempty(args.DualTree)
    
    if isempty(args.Image)
        error('You must supply either an image or dual-tree');
    else
        args.DualTree = dtwavexfm2(args.Image, max(args.Levels) + 1);
        [r c] = size(args.Image);
    end
else
    [r c] = size(args.DualTree{1}(:,:,1));
    r = 2*r;
    c = 2*c;
end

%Compute the ICP transformation of the image
[ilp] = mb_dual_tree_transform(args.DualTree);

%pre-allocate space for upscaled ICP map at each level
max_ilp_full = zeros(r,c,length(args.Levels));

%Compute upscaled ICP map at each level (could try and interpolate here)
for lev = 1:length(args.Levels)
    level = args.Levels(lev);
    temp = kron(max(ilp{level}, [], 3), ones(2^level));
    max_ilp_full(:,:,lev) = temp(1:r, 1:c);
end

%Compute maximum orientation across all levels
ilp_map = max(max_ilp_full, [], 3);