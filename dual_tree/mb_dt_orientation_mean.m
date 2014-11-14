function [orientation, orientation_levels] = mb_dt_orientation_mean(varargin)
%MB_DT_ORIENTATION_SUM *Insert a one line summary here*
%   [orientation,orientation_levels] = mb_dt_orientation_sum(varargin)
%
% MB_DT_ORIENTATION_SUM uses the U_PACKARGS interface function
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
    'Levels', 5);

%Check the user has supplied either an image or a dual_tree
if isempty(args.DualTree)
    
    if isempty(args.Image)
        error('You must supply either an image or dual-tree');
    else
        args.DualTree = dtwavexfm2(args.Image, args.Levels + 1);
    end
end

%pre-allocate space for the orientation mean in each level
orientation_levels = zeros(args.Levels, 2);

%Compute the ICP transformation of the image
[ilp icp] = mb_dual_tree_transform(args.DualTree);
clear ilp;

%Now compute the orientation sum in each level
for lev = 1:args.Levels
    
    %swap the orientation so everything lies between 0 and pi (otherwise
    %opposite directions (e.g. theta, theta+pi) would cancel each other out
    %NO!!! - This screws up the sum - instead report the mean from both
    %hemispheres
    
    max_icp = max(icp{lev}, [], 3);
    n_idx = imag(max_icp) >= 0;
    s_idx = imag(max_icp) < 0;
    
    %if we're using the whole image, just sum across all icp band in the
    %level
    if isempty(args.Subs)
        orientation_levels(lev, 1) = sum(max_icp(n_idx));
        orientation_levels(lev, 2) = sum(max_icp(s_idx));
    else
        %Down sample subs for this level
        subs = ceil(args.Subs / 2^lev);
        
        %convert to indices
        idx = sub2ind(size(icp{lev}(:,:,1)), subs(:,1), subs(:,2));
        
        %take mean of icp at indices across each band
        orientation_levels(lev, 1) = mean(max_icp(intersect(idx, n_idx)));
        orientation_levels(lev, 2) = mean(max_icp(intersect(idx, s_idx)));
    end
    
end

%Now combine sum adjusting for level size
orientation = mean(orientation_levels);