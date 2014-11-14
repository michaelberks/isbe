function [dual_tree] = mb_dual_tree_transform_i(ilp, icp)
%MB_DUAL_TREE_TRANSFORM_I *Insert a one line summary here*
%   [ilp icp] = mb_dual_tree_transform_i(dual_tree) - calculate the
%   inter-level (ILP) and inter-coefficient (ICP) products from a dual-tree
%   complex wavelet transform. The ILP will provide information on the type
%   of feature (eg. line vs edge) at any point, whilst the ICP provides the
%   orientation. See M Berks documentation
%
% Inputs:
%      ilp- *Insert description of input variable here*
%      icp-
%
%
% Outputs:
%      dual_tree- *Insert description of input variable here*
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 16-Oct-2008
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester

%Go through the tree from coarse to fine generating the dual-tree
%coefficients from the ILP and ICP coefficinets

num_levels = size(ilp,1);

dual_tree = cell(num_levels,1);
dual_tree{num_levels} = ilp{num_levels};

for level = num_levels-1:-1:1
    
    [Y X] = size(ilp{level}(:,:,1));
    dual_tree{level} = zeros(Y, X, 6);
    
    % Interpolate coarser level using complex interpolator.
    
    % Set up the expected phase shifts for each subband:
    % w = [-3 -1; -3 -3; -1 -3; 1 -3; 3 -3; 3 -1]*pi/2.15; 
    w = [-3 -1; -3 -3; -1 -3; 1 -3; 3 -3; 3 -1]*pi/2.15; % Nominally j * pi/2, but reduced a bit due to asymmetry of subband freq responses.
    p = [1 3]/4;  % Interpolation points
    
    for band = 1:6,
        % clipping extra space, for the moment
        temp = cpxinterp2(dual_tree{level+1}(:,:,band), p-0.5, w(band,:),'spline');        
        dt_int = temp(1:Y, 1:X);
        clear temp;
        
        % Double the angle of interpolated coefficients
        dt_int2 = dt_int .* dt_int ./ (abs(dt_int) + 1e-6);
        
        %Now we need check the phase of the ICP coefficient and negate the
        %real part of the ILP coefficient if the ICP angle is > pi
        ilp_phase = ilp{level}(:,:,band);
        idx = imag(icp{level}(:,:,band)) < 0;
        ilp_phase(idx) = complex(-real(ilp_phase(idx)), imag(ilp_phase(idx)));
        
        % Generate new dual-tree coefficients by summing the phases of the
        % phase-doubled, interpolated coarser level coefficients with the
        % ILP coefficients
        dt_phase = ilp_phase .* dt_int2;    
        
        %for sub-bands 4,5,6 rotate vector by pi/2
        if band >= 4
            dt_phase = complex(-imag(dt_phase), real(dt_phase));
        end
        
        %Take phase from angle sum and magnitude from ICP and copy into dual-tree structure        
        dual_tree{level}(:,:,band) = abs(icp{level}(:,:,band)).*exp(i*angle(dt_phase));
        
    end
    
end

%Finally invert the phase of 2nd and 5th sub-bands 
for level = 1:num_levels
    for band = [2 5]
        dual_tree{level}(:,:,band) = -dual_tree{level}(:,:,band);
    end
end