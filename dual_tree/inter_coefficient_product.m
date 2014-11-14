function [icp Wld] = inter_coefficient_product(dual_tree)
%INTERSCALE_COEFFICIENT_PRODUCT
%   [icp] = interscale_coefficient_product(dual_tree)
%
% Inputs:
%      dual_tree- pyramid structure of DT-CWT coefficients to transform
%
%
% Outputs:
%      icp- pyramid structure (matching the first n-1 levels of the dual-tree)
%           of inter-coefficient products
%       Wld- pyramid structure (matching the first n-1 levels of the dual-tree)
%           of coefficient phase differences
%
%
% Example:
%     [icp Wld] = inter_coefficient_product(dual_tree)  
%
% Notes:
%       For a full introduction on calculating the ILP and ICP see the
%       papers:
%           "Determining Multi-scale Image Feature Angles from Complex
%           Wavelet Phases" - Anderson, Kingsbury and Fauquer
%           "Coarse-level Object recognition using Inter-level Products of
%           Complex Wavelets" - Anderson, Kingsbury and Fauquer
%
% See also:
%       MB_DUAL_TREE_TRANSFORM and INTER_LEVEL_PRODUCT
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
    
    [Y X] = size(dual_tree{level}(:,:,1));
    if X < 2 || Y < 2 
        break;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % First do the ICP stuff
    %Create constant-phase complex values
    Wld{level} = zeros(Y, X, 6);
    
    x = 1:X - 1;
    y = 1:Y - 1;
    
    % Wl(x,y,1) = Wl(x,y,1)   .  Wl(x+1,y,1)*
    A = zeros(Y, X); B = zeros(Y, X);
    A(:,x) = dual_tree{level}(:,x,1) .* conj(dual_tree{level}(:,x+1,1));
    B(:,x+1) = A(:,x);
    Wld{level}(:,:,1) = max(A, B);
    
    % Wl(x,y,2) = Wl(x,y+1,2) .  Wl(x+1,y,2)*
    A = zeros(Y, X); B = zeros(Y, X);
    %A(y,x) = dual_tree{level}(y+1,x,2) .* conj(dual_tree{level}(y,x+1,2));
    A(y,x+1) = dual_tree{level}(y+1,x,2) .* conj(dual_tree{level}(y,x+1,2));
    B(y+1,x) = A(y,x+1);
    Wld{level}(:,:,2) = max(A, B);
    Wld{level}(1,1,2) = 1e-3; %Bit of a fudge here
    Wld{level}(Y,X,2) = 1e-3; %Bit of a fudge here
    
    % Wl(x,y,3) = Wl(x,y+1,3)   .  Wl(x,y,3)*
    A = zeros(Y, X); B = zeros(Y, X);
    A(y,:) = dual_tree{level}(y+1,:,3) .* conj(dual_tree{level}(y,:,3));
    B(y+1,:) = A(y,:);
    Wld{level}(:,:,3) = max(A, B);
    
    % Wl(x,y,4) = Wl(x,y,4)*  .  Wl(x,y+1,4)
    A = zeros(Y, X); B = zeros(Y, X);
    A(y,:) = conj(dual_tree{level}(y,:,4)) .* dual_tree{level}(y+1,:,4);
    B(y+1,:) = A(y,:);
    Wld{level}(:,:,4) = max(A, B);
    
    % Wl(x,y,5) = Wl(x,y,5)*  .  Wl(x+1,y+1,5)
    A = zeros(Y, X); B = zeros(Y, X);
    A(y,x) = conj(dual_tree{level}(y,x,5)) .* dual_tree{level}(y+1,x+1,5);
    B(y+1,x+1) = A(y,x);
    Wld{level}(:,:,5) = max(A, B);
    Wld{level}(1,X,5) = 1e-3; %Bit of a fudge here
    Wld{level}(Y,1,5) = 1e-3; %Bit of a fudge here
    
    % Wl(x,y,6) = Wl(x,y,1)*  .  Wl(x+1,y,6)
    A = zeros(Y, X); B = zeros(Y, X);
    A(:,x) =  conj(dual_tree{level}(:,x,6)) .*dual_tree{level}(:,x+1,6);
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
    
    %icp_mag{level} = abs(dual_tree{level}); %return DT magnitude
    icp_mag{level} = sqrt(abs(Wld{level})); %return ICP magnitude
    
    %Combine phase and magnitude into complex coefficients
    icp{level} = exp(i*icp_phase{level}) .* icp_mag{level};
end
    
    