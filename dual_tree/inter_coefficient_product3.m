function [icp Wld] = inter_coefficient_product3(dual_tree)
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
icp = cell(num_levels, 1);
Wld = cell(num_levels, 1);

for level = 1:num_levels
    
    if any(size(dual_tree{level}(:,:,1)) < 2)
        break;
    end
    
    %Create constant-phase complex values
    Wld{level} = zeros(size(dual_tree{level}));
    
    % Wl(x,y,1) = Wl(x,y,1)   .  Wl(x+1,y,1)*
    %          y    x                          y     x                               y  x+1                      
    Wld{level}(:,1:end-1,1) = dual_tree{level}(:,1:end-1,1) .* conj(dual_tree{level}(:,2:end,1));
    
    % Wl(x,y,2) = Wl(x,y+1,2) .  Wl(x+1,y,2)*
    %          y          x                          y+1     x                                   y      x+1       
    Wld{level}(1:end-1,1:end-1,2) = dual_tree{level}(2:end,1:end-1,2) .* conj(dual_tree{level}(1:end-1,2:end,2));
    
    % Wl(x,y,3) = Wl(x,y+1,3)   .  Wl(x,y,3)*
    %           y      x                        y+1  x                               y     x
    Wld{level}(1:end-1,:,3) = dual_tree{level}(2:end,:,3) .* conj(dual_tree{level}(1:end-1,:,3));
    
    % Wl(x,y,4) = Wl(x,y,4)  .  Wl(x,y+1,4)*
    %           y      x                              y     x                         y+1   x
    Wld{level}(1:end-1,:,4) = dual_tree{level}(1:end-1,:,4) .* conj(dual_tree{level}(2:end,:,4));
    
    % Wl(x,y,5) = Wl(x,y,5)  .  Wl(x+1,y+1,5)*
    %            y       x                                   y      x                             y+1   x+1
    Wld{level}(1:end-1,1:end-1,5) = dual_tree{level}(1:end-1,1:end-1,5) .* conj(dual_tree{level}(2:end,2:end,5));
    
    % Wl(x,y,6) = Wl(x,y,1)  .  Wl(x+1,y,6)*
    %          y   x                                y   x                             y  x+1
    Wld{level}(:,1:end-1,6) = dual_tree{level}(:,1:end-1,6) .* conj(dual_tree{level}(:,2:end,6));
    
    %Calculate inter-coefficient product phases
    icp_phase{level}(:,:,1) = (angle(Wld{level}(:,:,1)) / 4.49);
    icp_phase{level}(:,:,2) = (angle(Wld{level}(:,:,2)) / 8.98) + pi/4;
    icp_phase{level}(:,:,3) = (angle(Wld{level}(:,:,3)) / 4.49) + pi/2;
    icp_phase{level}(:,:,4) = (angle(Wld{level}(:,:,4)) / 4.49) + pi/2;
    icp_phase{level}(:,:,5) = (angle(Wld{level}(:,:,5)) / 8.98) + 3*pi/4;
    icp_phase{level}(:,:,6) = (angle(Wld{level}(:,:,6)) / 4.49);
    
    %Calculate inter-coefficient product magnitudes
    icp_mag{level} = abs(Wld{level});
    
    icp{level} = complex(cos(icp_phase{level}), sin(icp_phase{level})) .* icp_mag{level};
end
    
    