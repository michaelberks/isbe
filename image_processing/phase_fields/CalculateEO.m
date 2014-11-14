function [E_O]=CalculateEO(phi_x, gx, gy, alpha, lambda)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% 
% function CalculateEO
%
% syntax:
% 
% Output: 
%       output: the output image.
%            
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
if nargin < 4
    alpha = 0.0005;
    lambda = 0.5;
end

U_phi_x = lambda * ((0.25*phi_x.^4) - 0.5*phi_x.^2) + alpha * (phi_x - (phi_x.^3)/3);

E_O_phi_x = 0.5*(gx.^2 + gy.^2) + U_phi_x;
E_O = sum(E_O_phi_x(:));