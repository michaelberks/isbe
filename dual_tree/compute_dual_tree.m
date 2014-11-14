function [dual_tree] = compute_dual_tree(image_in, num_levels, use_nag)
%COMPUTE_DUAL_TREE compute dual tree decomposition of an image and either
%return the raw tree, or if use_nag is true, compute the interpolation
%structure of knot points using the NAG libraries
%   [dual_tree] = compute_dual_tree(image_in, num_levels, use_nag)
%
% Inputs:
%      image_in - image to decompose
%
%      num_levels - number of levels in the dual-tree
%
%      use_nag - flag specifying whether to return the raw tree or an
%      interpolation structure computed using the NAG libraries
%
%
% Outputs:
%      dual_tree - either the raw dual-tree or the interpolation structure
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 14-Sep-2010
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

if ~exist('use_nag','var'), use_nag = false; end

% Create DT-CWT of image
dt = dtwavexfm2b(image_in, num_levels);

%Try computing NAG interpolation - this may not work eg. if NAG libraries
%aren't installed
if use_nag
    try 
        [dual_tree.knot_mag dual_tree.knot_im dual_tree.knot_re dual_tree.dt_dims] = dt_interp_nag_in(dt);
    catch
        err = lasterror;
        display('Problem executing NAG interpolation. Will use matlab interpolation');
        display(err.message);
        dual_tree = dt;
    end
else
    dual_tree = dt;
end