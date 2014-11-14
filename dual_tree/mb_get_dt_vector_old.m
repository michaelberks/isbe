function [dt_coeffs] = mb_get_dt_vector_old(varargin)
%MB_GET_DT_VECTOR
%   [dt_coeffs] = mb_get_dt_vector(varargin)
%
% MB_GET_DT_VECTOR uses the U_PACKARGS interface function
% and so arguments to the function may be passed as name-value pairs
% or as a struct with fields with corresponding names. The field names
% are defined as:
%
% Mandatory Arguments:
%   - pts: a 2xN vector of (r,c) coordinates for the image level of a dual-tree
%          note these are at a factor of 2 higher resolution than level 1 of the
%
%   - dual_tree:
%
% % Optional Arguments:
%   - theta: {'weighted', maximum} choose which method you use to compute
%   the reference orientation at each level.
%       'weighted' - take sum of orientations across all scales weighted by
%          the icp amplitude at each scale
%       'maximum' - take the orientation of the scale with maximum icp amplitude
%
% Outputs:
%      dt_coeffs- vector of output coefficients computed for each location
%      from the input dual-tree: for each location we compute the follwoing
%      coefficients
%           - theta: reference orientation for the location either the
%           weighted sum of orientations across each level or the
%           orientation of the level with highest ampltiude
%           - A(i): for level i = 1:L, the amplitude of the original
%           dual-tree
%           -phi(i): for level i = 1:L, the orientation relative to theta
%           - delta(i): for level i = 1:L-1, the inter-scale phase
%           difference
%      As a result dt_coeffs is an Nx3*levels array - each corresponds to 
%      the coefficients for given input location. To make the indexing
%      neat, we put theta at the end of each row so we have
%      dt_coeffs(pt,:) = [A(1),phi(1),delta(1),...,A(L-1),phi(L-1),delta(L-1),A(L), phi(L), theta] 
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

% Unpack the arguments:
args = u_packargs(varargin, '0',...
            {... % Mandatory arguments
            'pts',...
            'dual_tree'...
            }, ... 
			'theta', 'weighted' ...
			);

% Get number of levels in tree
levels = length(args.dual_tree) - 1;

%Compute which sub-band has maximum amplitude at each location in each
%level of the DT
subband_idx = cell(levels, 1);
for lev = 1:levels
    [dummy subband_idx{lev}] = max(args.dual_tree{lev}, [], 3);
    clear dummy;
end    
    
% Compute the ICP and ILP transformations of the dual_tree
[ilp icp] = mb_dual_tree_transform(args.dual_tree);

% Convert theta method to logical input to avoid expensive string compares
% in main loop
if strcmpi(args.theta, 'weighted')
    weight_theta = true;
elseif strcmpi(args.theta, 'maximum');
    weight_theta = false;
else
    warning('Method for computing theta not recognised, using ''weighted method instead'); %#ok
    weight_theta = true;
end
    

% 
% obtain the vector path of parents coefficients through each subsequent
% level of the dual_tree
num_pts = size(args.pts, 1);

dt_coeffs = zeros(num_pts, 3*(levels-1) + 1);

%for now lets loop through each row of pts - in future we can make this
%more efficient
for pt = 1:num_pts
    
    r = args.pts(pt, 1); c = args.pts(pt, 2);
    
    icp_temp = zeros(levels, 1);
    for lev = 1:levels-1
        r_lev = ceil(r / 2^lev); %see resolution note above (r,c) are at image level
        c_lev = ceil(c / 2^lev);
        
        % get maximum ampitude
        dt_coeffs(pt, 3*lev - 2) = abs(args.dual_tree{lev}(r_lev, c_lev, subband_idx{lev}(r_lev, c_lev)));
        
        % get ICP coefficient - we need to make these relative;
        icp_temp(lev) = icp{lev}(r_lev, c_lev, subband_idx{lev}(r_lev, c_lev));
        
        % get interscale phase
        dt_coeffs(pt, 3*lev) = angle(ilp{lev}(r_lev, c_lev, subband_idx{lev}(r_lev, c_lev)));
       
    end
    
    %Now compute reference orientation theta
    if weight_theta
        theta = sum(icp_temp);
    else
        theta = max(icp_temp); %max acts on abs, then further by phase if necessary
    end
    
    %make phi(i) relative to theta, and make sure we're in space [-pi,pi)
    %by using multiplication of complex conjugates then taking angle
    dt_coeffs(pt, 2:3:end-1) = angle(icp_temp .* conj(theta));
    
    % Put theta at the end of the row of coefficients
    dt_coeffs(pt, end) = theta;
    
end
    
