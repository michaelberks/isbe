function [icp_full] = mb_full_image_icp(image_in, min_level, max_level)
%DT_PHASE_CONGRUENCY *Insert a one line summary here*
%   [phase_map] = dt_phase_congruency(image_in,levels)
%
% Inputs:
%      image_in- *Insert description of input variable here*
%
%      levels- *Insert description of input variable here*
%
%
% Outputs:
%      phase_map- *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 07-May-2009
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester

[Y X] = size(image_in);

%build dual-tree
dual_tree = dtwavexfm2(image_in, max_level+1);

%compute ICP
[ilp icp] = mb_dual_tree_transform(dual_tree); clear ilp;

%convert to the full tree and take magnitudes (we don't need the full tree
%phase)
dual_tree_full_abs = abs(dt_to_full_image(dual_tree, min_level:max_level));
clear dual_tree;

%pre-allocate full icp structure
%icp_full = zeros(Y, X, size(dual_tree_full_abs,4));
icp_full = zeros(size(dual_tree_full_abs));

%Upscale the icp coefficients (don't interpolate) and compute the maximum
%band from which to sample at each pixel
for lev = max_level-min_level+1:-1:1
    level = min_level + lev - 1;
    
    %Get indices to maximal bands
    %[dummy max_dt_bands] = max(dual_tree_full_abs(:,:,:,lev), [], 3); clear dummy
    
    %Pre-allocate icp output for this level
    %icp_level = zeros(Y, X);
    icp_level = zeros(Y, X, 6);
    
    for band = 1:6
        
        %Workout pixels maximal in this band/level
        %band_idx = max_dt_bands == band;
        
        %Upscale ICP coefficients and compute phase
        icp_band_phase = kron(angle(icp{level}(:,:,band)), ones(2^level));
        
        %Get magnitudes from interpolated DT-CWT magnitudes
        dt_band_abs = dual_tree_full_abs(:,:,band,lev);
        
        %Combine ICP phase and DT magnitude
%         icp_level(band_idx) = ...
%             dt_band_abs(band_idx) .* exp(i*icp_band_phase(band_idx));
        icp_level(:,:,band) = dt_band_abs .* exp(i*icp_band_phase);
    end
    %icp_full(:,:,lev) = icp_level;
    icp_full(:,:,:,lev) = icp_level;
end
%icp_full(imag(icp_full) < 0) = -icp_full(imag(icp_full) < 0);

% %convert to the full tree
% icp = dt_to_full_image(dual_tree, min_level:max_level);
% clear dual_tree;
% 
% %Do the ILP stuff - i.e. for each level from coarse to fine subtract the
% %phase doubled coarser level
% for lev = max_level-min_level+1:-1:1
%     level = min_level + lev - 1;
%     h = 1;%2^(level - 1);
%     
%     x = h+1:X-h;
%     y = h+1:X-h;
%     
%     %do icp stuff
%     % Wl(x,y,1) = Wl(x-1,y,1)   .  Wl(x+1,y,1)*
%     icp(y,x,1,lev) = icp(y, x-h, 1,lev) .* conj(icp(y, x+h, 1,lev));
%     
%     % Wl(x,y,2) = Wl(x-1,y+1,2) .  Wl(x+1,y-1,2)*
%     icp(y,x,2,lev) = icp(y+h, x-h, 2,lev) .* conj(icp(y-h, x+h, 2,lev));
%     
%     % Wl(x,y,3) = Wl(x,y+1,3)   .  Wl(x,y-1,3)*
%     icp(y,x,3,lev) = icp(y+h, x, 3,lev) .* conj(icp(y-h, x, 3,lev));
%     
%     % Wl(x,y,4) = Wl(x,y-1,4)*  .  Wl(x,y+1,4)
%     icp(y,x,4,lev) = conj(icp(y-h, x, 4,lev)) .* icp(y+h, x, 4,lev);
%     
%     % Wl(x,y,5) = Wl(x-1,y-1,5)*  .  Wl(x+1,y+1,5)
%     icp(y,x,5,lev) = conj(icp(y-h, x-h, 5,lev)) .* icp(y+h, x+h, 5,lev);
%     
%     % Wl(x,y,6) = Wl(x-1,y,1)*  .  Wl(x+1,y,6)
%     icp(y,x,6,lev) = conj(icp(y, x-h, 6,lev)) .* icp(y, x+h, 6,lev);
%     
%     icp(:,:,1,lev) = sqrt(abs(icp(:,:,1,lev))).*exp(i*((2^level)*angle(icp(:,:,1,lev)) / 4.49));
%     icp(:,:,2,lev) = sqrt(abs(icp(:,:,2,lev))).*exp(i*((2^level)*angle(icp(:,:,2,lev)) / 8.98) + pi/4);
%     icp(:,:,3,lev) = sqrt(abs(icp(:,:,3,lev))).*exp(i*((2^level)*angle(icp(:,:,3,lev)) / 4.49) + pi/2);
%     icp(:,:,4,lev) = sqrt(abs(icp(:,:,4,lev))).*exp(i*((2^level)*angle(icp(:,:,4,lev)) / 4.49) + pi/2);
%     icp(:,:,5,lev) = sqrt(abs(icp(:,:,5,lev))).*exp(i*((2^level)*angle(icp(:,:,5,lev)) / 8.98) + 3*pi/4);
%     icp(:,:,6,lev) = sqrt(abs(icp(:,:,6,lev))).*exp(i*((2^level)*angle(icp(:,:,6,lev)) / 4.49));
%     
%     icp([1:h Y-h+1:Y],:,:,lev) = 1e-4;
%     icp(:,[1:h X-h+1:X],:,lev) = 1e-4;
%      
% end
% 
% icp(imag(icp) < 0) = -icp(imag(icp) < 0);
        
    
    
