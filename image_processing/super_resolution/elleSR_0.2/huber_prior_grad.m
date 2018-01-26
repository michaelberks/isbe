function G = huber_prior_grad(Dx, alpha, ni, nj)
%ELLE_EVAL_HUBER_GRAD
%   [] = elle_eval_huber_grad()
%
% Inputs:
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 28-Feb-2017
% Author: Michael Berks 
% Email : michael.berks@postgrad.man.ac.uk 
% Phone : +44 (0)161 275 5285 
% Copyright: (C) University of Manchester
%Dx is already |Dx|
rho_grad = hubergradx(Dx, alpha);

% huber_prior_gradx(G, X, biv, bih, nu) */
% adds the gradient of p(x) wrt x to the current values of the elements of G. */
G =zeros(ni*nj,1);
r2 = sqrt(2);

% loop to go back from vectors of size (4n)*(n) to vectors of size n*n. */
for i = 1:ni %(i=0; i<biv; i++) {
    for j = 1:nj %(j=0; j<bih; j++) {
        c = (j-1)*ni+i;
        tsum = 0;
        if i>1
           tsum = tsum + rho_grad(c-1, 1);
        end
        if j>1 
           tsum = tsum + rho_grad(c-ni, 2);
        end
        if i>1 && j>1
           tsum = tsum + rho_grad(c-ni-1, 3)/r2;
        end
        if i<ni && j>1 
           tsum = tsum + rho_grad(c-ni+1, 4)/r2;
        end
        if i<ni
           tsum = tsum - rho_grad(c, 1);
        end
        if j<nj
           tsum = tsum - rho_grad(c, 2);
        end
        if j<nj && i<ni
           tsum = tsum - rho_grad(c, 3)/r2;
        end
        if j<nj && i>1 
           tsum = tsum - rho_grad(c, 4)/r2;
        end

        G(c) = tsum;
    end
end


%----------------------------------------------------------------------
function [grad] = hubergradx(x, alpha)

    grad = 2*min(x, alpha);
    