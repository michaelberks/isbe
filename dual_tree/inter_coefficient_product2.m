function [icp] = inter_coefficient_product2(dual_tree)
%INTERSCALE_COEFFICIENT_PRODUCT *Insert a one line summary here*
%   [icp] = interscale_coefficient_product(dual_tree)
%
% Inputs:
%      dual_tree- *Insert description of input variable here*
%
%
% Outputs:
%      icp- *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 15-Oct-2008
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 1241 
% Copyright: (C) University of Manchester 

num_levels = length(dual_tree) - 1;
icp_mag = cell(num_levels, 1);
icp_phase = cell(num_levels, 1);
icp_phase2 = cell(num_levels, 1);
icp = cell(num_levels, 1);

for level = 1:num_levels
    
    [Y X] = size(dual_tree{level}(:,:,1));
    
    if X < 2 || Y < 2 
        break;
    end
    
    dti = zeros(2*Y + 1, 2*X + 1, 6);
    %dti = zeros(Y, X, 6);
    %w = [-3 -1; -3 -3; -1 -3; 1 -3; 3 -3; 3 -1]*pi/2.15;
    
    for ori = 1:6       
       dti(:,:,ori) = ...
           interp2(angle(dual_tree{level}(:,:,ori)), (1:2*X+1)/2, (1:2*Y+1)'/2, 'spline');
           %angle(dual_tree{level}(:,:,ori));
           %complex_interp2(dual_tree{level}(:,:,ori), (1:2*X+1)/2, (1:2*Y+1)/2, w(ori,:), 'spline');
    end

    x = 2:2:2*X;
    y = 2:2:2*Y;
    %x = 2:X-1;
    %y = 2:Y-1;
    
    %Create constant-phase complex values
    icp_phase{level} = zeros(size(dual_tree{level}));
    icp_phase2{level} = zeros(size(dual_tree{level}));
    icp_mag{level} = zeros(size(dual_tree{level}));
    
    % Wl(x,y,1) = Wl(x,y,1)   .  Wl(x+1,y,1)*                     
    %icp_phase{level}(y,x,1) = (dti(y,x,1) - dti(y,x+1,1)); 
    icp_phase2{level}(:,:,1) = (dti(y,x-1,1) - dti(y,x+1,1)); 
    
    % Wl(x,y,2) = Wl(x,y+1,2) .  Wl(x+1,y,2)*
    %icp_phase{level}(y,x,2) = (dti(y+1,x,2) - dti(y,x+1,2));% / 8.98); + pi/4;  
    icp_phase2{level}(:,:,2) = (dti(y+1,x-1,2) - dti(y-1,x+1,2));
    
    % Wl(x,y,3) = Wl(x,y+1,3)   .  Wl(x,y,3)*
    %icp_phase{level}(y,x,3) = (dti(y+1,x,3) - dti(y,x,3));% / 4.49); + pi/2;
    icp_phase2{level}(:,:,3) = (dti(y+1,x,3) - dti(y-1,x,3));
    
    % Wl(x,y,4) = Wl(x,y,4)*  .  Wl(x,y+1,4)
    %icp_phase{level}(y,x,4) = (dti(y+1,x,4) - dti(y,x,4));% / 4.49);% + pi/2;
    icp_phase2{level}(:,:,4) = (dti(y+1,x,4) - dti(y-1,x,4));
    
    % Wl(x,y,5) = Wl(x,y,5)*  .  Wl(x+1,y+1,5)
    %icp_phase{level}(y,x,5) = (dti(y+1,x+1,5) - dti(y,x,5));% / 8.98); + 3*pi/4;
    icp_phase2{level}(:,:,5) = (dti(y+1,x+1,5) - dti(y-1,x-1,5));
    
    % Wl(x,y,6) = Wl(x,y,1)*  .  Wl(x+1,y,6)
    %icp_phase{level}(y,x,6) = (dti(y,x+1,6) - dti(y,x,6));% / 4.49;
    icp_phase2{level}(:,:,6) = (dti(y,x+1,6) - dti(y,x-1,6));
    
    %Calculate inter-coefficient product phases
%     icp_phase{level}(:,:,1) = (pi_rotate(icp_phase{level}(:,:,1)) / 4.49);
%     icp_phase{level}(:,:,2) = (pi_rotate(icp_phase{level}(:,:,2)) / 8.98) + pi/4;
%     icp_phase{level}(:,:,3) = (pi_rotate(icp_phase{level}(:,:,3)) / 4.49) + pi/2;
%     icp_phase{level}(:,:,4) = (pi_rotate(icp_phase{level}(:,:,4)) / 4.49) + pi/2;
%     icp_phase{level}(:,:,5) = (pi_rotate(icp_phase{level}(:,:,5)) / 8.98) + 3*pi/4;
%     icp_phase{level}(:,:,6) = (pi_rotate(icp_phase{level}(:,:,6)) / 4.49);
    
    icp_phase{level}(:,:,1) = (pi_rotate(icp_phase2{level}(:,:,1)) / 4.49);
    icp_phase{level}(:,:,2) = (pi_rotate(icp_phase2{level}(:,:,2)) / 8.98) + pi/4;
    icp_phase{level}(:,:,3) = (pi_rotate(icp_phase2{level}(:,:,3)) / 4.49) + pi/2;
    icp_phase{level}(:,:,4) = (pi_rotate(icp_phase2{level}(:,:,4)) / 4.49) + pi/2;
    icp_phase{level}(:,:,5) = (pi_rotate(icp_phase2{level}(:,:,5)) / 8.98) + 3*pi/4;
    icp_phase{level}(:,:,6) = (pi_rotate(icp_phase2{level}(:,:,6)) / 4.49);
    
    %Calculate inter-coefficient product magnitudes
    %icp_mag{level} = abs(Wld{level});
    %icp_mag{level}(y,x,:) = abs(dual_tree{level}(y,x,:));
    icp_mag{level} = abs(dual_tree{level});
    icp{level} = exp(i*icp_phase{level}) .* icp_mag{level};
end
    
    