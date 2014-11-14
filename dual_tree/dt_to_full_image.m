function [full_tree] = dt_to_full_image(dual_tree, levels, interpmethod)

%DT_TO_FULL_IMAGE interpolate each level of a DT-CWT to the full pixel grid
%   [full_tree] = dt_to_full_image(dual_tree)
%
% Inputs:
%      dual_tree- DT-CWT cell array
%
%      levels - 1D array specifying the levels of the dual-tree to include
%      in the full tree (to save memory it is useful not include the 1st
%      level for example)
%
%      interpmethod - string description of interpolation method passed to
%      Matlab's interp2 function
%
% Outputs:
%      full_tree - (m, n, 6, num_levels)-array such that (:,:,b,l) is the
%      array of interpolated dual-tree coefficients for b-th oriented
%      sub-band in the levels(l)-th level of the original dual tree
%
%
% Example:
%
% Notes:
%
% See also: COMPLEX_INTERP_FULL_IMAGE, CPXINTERP2
%
% Created: 05-May-2009
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester

%set default for interpmethod
if nargin < 3
    interpmethod = 'cubic';
end

if nargin < 2
    levels = 1:(size(dual_tree, 1) - 1);
end

%Get number of levels to put in full tree
num_levels = length(levels);

%Get size of full pixel grid (assume double the size of the top dual-tree
%level)
[m n] = size(dual_tree{1}(:,:,1));
M = 2*m;
N = 2*n;

%Set central frequencies of bands for interp method
% w = [-3 -1; -3 -3; -1 -3; 1 -3; 3 -3; 3 -1]*pi/2.15;
w1 = 16*0.0644;
w2 = 16*0.2859;
w3 = 16*0.2238;
w = [-w2 -w1; -w3 -w3; -w1 -w2; w1 -w2; w3 -w3; w2 -w1];
%w = [-w2 -w1; -w2 -w2; -w1 -w2; w1 -w2; w2 -w2; w2 -w1];

%Pre-allocate space for full-tree - can be an array since each level has
%same number of coeffs
full_tree = zeros(M, N, 6, num_levels);

%for each sub-band, interpolate up the coefficients to the full size
for l = 1:num_levels
    for band = 1:6
        full_tree(:,:,band,l) =...
            complex_interp_full_image(dual_tree{levels(l)}(:,:,band),...
            M, N, levels(l), w(band,:), interpmethod);
    end
end

%finally, apply the phase correction to each of the sub-bands - multiply
%bands by - note this still needs testing (it affects the ILP calculations)
% NK suggests {1,-1,1,j,-j,j}, but I think it should be {1,-1,1,-j,j,-j} -
% after checking I agree with NK!!!!
full_tree(:,:,[1 3],:) =  i*full_tree(:,:,[1 3],:);
full_tree(:,:,[4 6],:) =   -full_tree(:,:,[4 6],:);
full_tree(:,:,2,:)     = -i*full_tree(:,:,2,:);

full_tree(:,:,[1 3 4 6],1) =   -full_tree(:,:,[1 3 4 6],1);