function [dt_coeffs] = mb_get_dt_vector(pts, dual_tree, num_levels)
%MB_GET_DT_VECTOR
%   [dt_coeffs] = mb_get_dt_vector(varargin)
%
% Inputs:
%   - pts: a 2xN vector of (r,c) coordinates for the image level of a dual-tree
%          note these are at a factor of 2 higher resolution than level 1 of the
%
%   - dual_tree:
%
%
% Outputs:
%      dt_coeffs- vector of output coefficients computed for each location
%      from the input dual-tree: for each location we compute the follwoing
%      coefficient
%      
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 31-Oct-2008
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester 

% Get number of levels in tree
if nargin < 3
    num_levels = length(dual_tree) - 2;
end

% Compute the ICP and ILP transformations of the dual_tree
[ilp] = mb_dual_tree_transform(dual_tree);

%Compute which sub-band has maximum amplitude at each location in each
%level of the DT
subband_idx = cell(num_levels, 1);
for lev = 1:num_levels
    %[dummy subband_idx{lev}] = max(dual_tree{lev}, [], 3);
    [dummy subband_idx{lev}] = max(ilp{lev}, [], 3);
    clear dummy;
end    
    
% 
% obtain the vector path of parents coefficients through each subsequent
% level of the dual_tree
num_pts = size(pts, 1);

dt_coeffs = zeros(num_pts, 2*num_levels);    

%for now lets loop through each row of pts - in future we can make this
%more efficient
for pt = 1:num_pts
    
    r = pts(pt, 1); c = pts(pt, 2);
    
    for lev = 1:num_levels
        r_lev = ceil(r / 2^(lev-1)); %see resolution note above (r,c) are at image level
        c_lev = ceil(c / 2^(lev-1));
        
        % get maximum amplitude
        %dt_mag = abs(dual_tree{lev}(r_lev, c_lev, subband_idx{lev}(r_lev, c_lev)));
        dt_mag = abs(ilp{lev}(r_lev, c_lev, subband_idx{lev}(r_lev, c_lev))); %take ILP not DT magnitude
        % get interscale phase
        ilp_phase = angle(ilp{lev}(r_lev, c_lev, subband_idx{lev}(r_lev, c_lev)));
        
        dt_coeffs(pt, 2*lev - 1) = dt_mag;% * cos(ilp_phase);
        dt_coeffs(pt, 2*lev) = ilp_phase;%dt_mag * sin(ilp_phase);    
    end    
end   
