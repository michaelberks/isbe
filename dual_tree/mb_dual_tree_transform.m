function [ilp icp] = mb_dual_tree_transform(dual_tree, finest_level, coarsest_level)
%MB_DUAL_TREE_TRANSFORM
%   [ilp icp] = mb_dual_tree_transform(dual_tree) - calculate the
%   inter-level (ILP) and inter-coefficient (ICP) products from a dual-tree
%   complex wavelet transform. The ILP will provide information on the type
%   of feature (eg. line vs edge) at any point, whilst the ICP provides the
%   orientation. See notes below
%
% Inputs:
%      dual_tree- pyramid structure of DT-CWT coefficients to transform
%
%
% Outputs:
%      ilp- The pyramid structure of inter-level products
%      icp- The pyramid structure of inter-coefficients products
%
% Example:
%       [ilp icp] = mb_dual_tree_transform(dual_tree);
%
% Notes:
%       For a full introduction on calculating the ILP and ICP see the
%       papers:
%           "Determining Multi-scale Image Feature Angles from Complex
%           Wavelet Phases" - Anderson, Kingsbury and Fauquer
%           "Coarse-level Object recognition using Inter-level Products of
%           Complex Wavelets" - Anderson, Kingsbury and Fauquer
%
%       This code also implements our own modifications to both products
%       which will be documented in M Berks' thesis
%
% See also:
%       INTER_COEFFICIENT_PRODUCT, INTER_LEVEL_PRODUCT
%
% Created: 16-Oct-2008
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester
num_levels = size(dual_tree,1)-1;

if nargin < 2
    finest_level = [];
end
if nargin < 3
    coarsest_level = [];
end
if isempty(finest_level)
    finest_level = 1;
end
if isempty(coarsest_level)
    coarsest_level = num_levels-1;
end

%Change phase of 2nd and 5th subbands
for level = finest_level:coarsest_level+1
    dual_tree{level}(:,:,2) =   -dual_tree{level}(:,:,2);
	dual_tree{level}(:,:,4) =  i*dual_tree{level}(:,:,4);
	dual_tree{level}(:,:,5) = -i*dual_tree{level}(:,:,5);
	dual_tree{level}(:,:,6) =  i*dual_tree{level}(:,:,6);
end

%Pre-allocate ILP and ICP structures, storing the coarsest dual-tree
%wavelet level (i.e. the level above the scaling coefficients) in the ILP
%structure
ilp = cell(coarsest_level+1, 1);
ilp{coarsest_level+1} = dual_tree{coarsest_level+1};
icp = cell(coarsest_level+1, 1);

%Pre-allocate structures needed to compute ICP
icp_mag = cell(coarsest_level+1, 1);
icp_phase = cell(coarsest_level+1, 1);

%Corresponds to the W_l_delta matrix in the Anderson paper
Wld = cell(coarsest_level+1, 1);

%Work through the levels of the dual-tree (we go coarse-to-fine but this is
%an arbitrary choice)
for level = coarsest_level:-1:finest_level
    
    [Y X] = size(dual_tree{level}(:,:,1));
    if X < 2 || Y < 2 
        break;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % First do the ICP stuff
    %Create constant-phase complex values
    
    %Preallocate space for wavelet phase differences
    % We want to avoid zero value coefficients which we'll have in the top
    % left and bottom right pixels of bands 2/5 - so stick a small value in...
    Wld{level} = ones(Y, X, 6)*1e-4;
    
    x = 1:X - 1;
    y = 1:Y - 1;
    
    % For each sub-band in turn compute the phase difference between a
    % coefficient and its two neighbours, by multplying the coefficient
    % with the complex conjugate of its neighbours. We then take the
    % maximum of these two products and store in Wld. The directions of the
    % neighbours in each subband are: 1,6 - horizontal; 3,4 = vertical, 2 =
    % diagonal at 45 degrees; 5 = diagonal at 135 degrees
    
    % Wl(x,y,1) = Wl(x,y,1)   .  Wl(x+1,y,1)*
    A = zeros(Y, X); B = zeros(Y, X);
    A(:,x) = dual_tree{level}(:,x,1) .* conj(dual_tree{level}(:,x+1,1));
    B(:,x+1) = A(:,x);
    Wld{level}(:,:,1) = max(A, B);
    
    % Wl(x,y,2) = Wl(x,y+1,2) .  Wl(x+1,y,2)*
    A = zeros(Y, X); B = zeros(Y, X);
    A(y,x+1) = dual_tree{level}(y+1,x,2) .* conj(dual_tree{level}(y,x+1,2));
    B(y+1,x) = A(y,x+1); 
    Wld{level}(:,:,2) = max(A, B);
    
    % Wl(x,y,3) = Wl(x,y+1,3)   .  Wl(x,y,3)*
    A = zeros(Y, X); B = zeros(Y, X);
    A(y,:) = dual_tree{level}(y+1,:,3) .* conj(dual_tree{level}(y,:,3));
    B(y+1,:) = A(y,:);
    Wld{level}(:,:,3) = max(A, B);
    
    % Wl(x,y,4) = Wl(x,y,4)*  .  Wl(x,y+1,4)
    A = zeros(Y, X); B = zeros(Y, X);
    A(y,:) = dual_tree{level}(y+1,:,4) .* conj(dual_tree{level}(y,:,4));
    B(y+1,:) = A(y,:);
    Wld{level}(:,:,4) = max(A, B);
    
    % Wl(x,y,5) = Wl(x,y,5)*  .  Wl(x+1,y+1,5)
    A = zeros(Y, X); B = zeros(Y, X);
    A(y,x) = dual_tree{level}(y+1,x+1,5) .* conj(dual_tree{level}(y,x,5));
    B(y+1,x+1) = A(y,x);
    Wld{level}(:,:,5) = max(A, B);
    
    % Wl(x,y,6) = Wl(x,y,1)*  .  Wl(x+1,y,6)
    A = zeros(Y, X); B = zeros(Y, X);
    A(:,x) =  dual_tree{level}(:,x+1,6) .*conj(dual_tree{level}(:,x,6));
    B(:,x+1) = A(:,x);
    Wld{level}(:,:,6) = max(A, B);
    
    %Calculate inter-coefficient product phases, put in range [0, pi)
    icp_phase{level}(:,:,1) = mod((angle(Wld{level}(:,:,1)) / 4.49), pi);
    icp_phase{level}(:,:,2) = mod((angle(Wld{level}(:,:,2)) / 8.98) + pi/4, pi);
    icp_phase{level}(:,:,3) = mod((angle(Wld{level}(:,:,3)) / 4.49) + pi/2, pi);
    icp_phase{level}(:,:,4) = mod((angle(Wld{level}(:,:,4)) / 4.49) + pi/2, pi);
    icp_phase{level}(:,:,5) = mod((angle(Wld{level}(:,:,5)) / 8.98) + 3*pi/4, pi);
    icp_phase{level}(:,:,6) = mod((angle(Wld{level}(:,:,6)) / 4.49), pi);
    
     %Calculate inter-coefficient product magnitudes
    
    % I think it's better to maintain the magnitude from the original
    % dual-tree coefficient (although one could argue since we have this
    % already we may as well return the product coefficient from this
    % function - in which case (un)comment the statements below)
    
    icp_mag{level} = abs(dual_tree{level});
    %icp_mag{level} = sqrt(abs(Wld{level}));
    
    icp{level} = exp(i*icp_phase{level}) .* icp_mag{level};
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Now do the ILP stuff updating the ICP phase as we go
    ilp{level} = zeros(Y, X, 6);
    
    % Interpolate coarser level using complex interpolator.
    
    % Set up the expected phase shifts for each subband:
    w = [-3 -1; -3 -3; -1 -3; 1 -3; 3 -3; 3 -1]*pi/2.15; % Nominally j * pi/2, but reduced a bit due to asymmetry of subband freq responses.
    p = [1 3]/4;  % Interpolation points
    
    for band = 1:6,
        % clipping extra space, for the moment
        temp = cpxinterp2(dual_tree{level+1}(:,:,band), p-0.5, w(band,:),'spline');        
        dt_int = temp(1:Y, 1:X);
        clear temp;
        
        % Double the angle of interpolated coefficients
        %dt_int2 = dt_int .* dt_int ./ (abs(dt_int) + 1e-6);
        dt_int2 = abs(dt_int) .* exp(i*2*angle(dt_int));
        
        % Generate vectors at difference angles of current level coefficients and
        % phase-doubled, interpolated coarser level coefficients
        angle_diff = dual_tree{level}(:,:,band) .* conj(dt_int2);
        
        %Now map the 2nd and 3rd complex quadrants onto the 1st and 4th by
        %flipping coefficients with negative real part
        idx = real(angle_diff) < 0; %locations to swap
        angle_diff(idx) =...
            complex(-real(angle_diff(idx)), imag(angle_diff(idx)));
        
        %Copy into output ILP structure
        %ilp{level}(:,:,band) = exp(i*angle(angle_diff)) .* abs(dual_tree{level}(:,:,band));
        ilp{level}(:,:,band) = angle_diff;
        
        %subtract pi from ICP phase for swapped ILP coefficients by taking
        %the negative of the complex coefficients
        if 1
            % First make sure icp is defined on 0, pi 
            icp_band = icp{level}(:,:,band);
            icp_band(imag(icp_band) < 0) = -icp_band(imag(icp_band) < 0);
            % Now swap as required
            icp_band(idx) = -icp_band(idx);
            icp{level}(:,:,band) = icp_band;
        end
        
    end
    
end