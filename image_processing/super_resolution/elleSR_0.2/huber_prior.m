function [sum_rho] = huber_prior(Dx, alpha)
% eout = huber_prior(x,alpha,bet,biv,bih); */

%Dx = abs(Dx);
mask = Dx < alpha;
sum_rho = sum(Dx(mask).^2) + sum(2*Dx(~mask)*alpha - alpha^2); 


%-------------------------------------------------------------------------
% function rho = huber(x, alpha)
% if (x<0) 
%     x=-x;
% end
% if (x<alpha) 
%     rho = x.*x; 
% else
%     rho = (2*x*alpha - alpha*alpha);
% end

