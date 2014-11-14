function [ilp] = inter_level_product(dual_tree)
%INTER_LEVEL_PRODUCT
%   [ilp] = inter_level_product(dual_tree)
%
% Inputs:
%      dual_tree- pyramid structure of DT-CWT coefficients to transform
%
%
% Outputs:
%      ilp- pyramid structure (matching the first n-1 levels of the dual-tree)
%           of inter-level products
%
%
% Example:
%       [ilp] = inter_level_product(dual_tree)
%
% Notes:
%       For a full introduction on calculating the ILP and ICP see the
%       papers:
%           "Determining Multi-scale Image Feature Angles from Complex
%           Wavelet Phases" - Anderson, Kingsbury and Fauquer
%           "Coarse-level Object recognition using Inter-level Products of
%           Complex Wavelets" - Anderson, Kingsbury and Fauquer
% See also:
%       MB_DUAL_TREE_TRANSFORM and INTER_COEFFICIENT_PRODUCT
% Created: 16-Oct-2008
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester

%Get number of levels in pyramid
num_levels = length(dual_tree) - 1;

%Pre-allocate ILP output structure
ilp = cell(num_levels-1, 1);

%Compute ILP coefficients for each level of the tree
for level = 1:num_levels-1
    
    
    [yf xf] = size(dual_tree{level}(:,:,1));
    
    coarse_level_int = zeros(yf, xf, 6);
    
    %Set up inputs for interpolation - Nick Kingsbury's method
    %NK - Nominally j * pi/2, but reduced a bit due to asymmetry of subband freq responses.
    w = [-3 -1; -3 -3; -1 -3; 1 -3; 3 -3; 3 -1]*pi/2.15; 
    p = [-.25 .25];  % Interpolation points
    
    for ori = 1:6  
        
        %Interpolate coarser level - NK's method     
        temp = cpxinterp2(dual_tree{level+1}(:,:,ori), p, w(ori,:),'spline');
        band_coarse = temp(1:yf, 1:xf);
        
        %Double the phase of band_coarse
        coarse_level_int(:,:,ori) = abs(band_coarse).*exp(2*angle(band_coarse)*i);   
        
    end
    %Multiply fine coefficients by the conjugates of the phase-doubled
    %coarse coefficients
    ilp{level} = dual_tree{level} .* conj(coarse_level_int);
end
