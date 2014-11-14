function [missing_mask] = wrong_mask_recovery(correct_mask, wrong_mask)
%WRONG_MASK_RECOVERY what to do if you predicted an image using the wrong
%mask!!
%   [missing_mask] = wrong_mask_recovery(correct_mask, wrong_mask, max_size)
%
% Inputs:
%      correct_mask - the mask you should have used
%
%      wrong_mask - the mask you did use
%
%
% Outputs:
%      missing_mask - a mask set to 1 at all pixels you didn't predict
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 27-Jan-2010
% Author: Michael Berks
% Email : michael.berks@postgrad.man.ac.uk
% Phone : +44 (0)161 275 1241
% Copyright: (C) University of Manchester

%Get size of correct mask and reshape the wrong mask to match (as predict_image_set would do)
[ROW COL] = size(correct_mask);
wrong_mask = imresize(wrong_mask, [ROW COL]);

missing_mask = correct_mask & ~wrong_mask;
