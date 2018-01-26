function G = elle_eval_huber_grad(X, alpha, beta, biv, bih)
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
if nargin < 3
    error('Wrong number of inputs (Should have: X, alpha, beta, biv, bih)');
end

% huber_prior_gradx(G, X, biv, bih, nu) */
% adds the gradient of p(x) wrt x to the current values of the elements of G. */
G =zeros(1,biv*bih);
   
% huber_prior_gradx(G, X, biv, bih, nu) */
% adds the gradient of p(x) wrt x to the current values of the elements of G. */
Y2 = zeros(biv*bih*4,1); % This is initialised to all zeros. */
oset1 = biv*bih;
oset2 = 2*oset1;
oset3 = oset1+oset2;
r2 = sqrt(2);

for i = 1:biv %(i=0; i<biv; i++) {
    for j = 1:bih %(j=0; j<bih; j++) {
        c = (j-1)*biv+i;
        if i<biv 
            Y2(c) = hubergradx(-X(c)+X(c+1), alpha);
            if j<bih 
                Y2(c+oset2) = hubergradx((-X(c)+X(c+biv+1))/r2, alpha); 
            end
        end
        if j<bih
            Y2(c+oset1) = hubergradx(-X(c)+X(c+biv), alpha);
            if i >= 2  
                Y2(c+oset3) = hubergradx((-X(c)+X(c+biv-1))/r2, alpha); 
            end
        end
    end
end

% loop to go back from vectors of size (4n)*(n) to vectors of size n*n. */
for i = 1:biv %(i=0; i<biv; i++) {
    for j = 1:bih %(j=0; j<bih; j++) {
        c = (j-1)*biv+i;
        tsum = 0;
        if i>1
           tsum = tsum + Y2(c-1);
        end
        if j>1 
           tsum = tsum + Y2(c-biv+oset1);
        end
        if i>1 && j>1
           tsum = tsum + Y2(c-biv-1+oset2)/r2;
        end
        if i<biv && j>1 
           tsum = tsum + Y2(c-biv+1+oset3)/r2;
        end
        if i<biv
           tsum = tsum - Y2(c);
        end
        if j<bih
           tsum = tsum - Y2(c+oset1);
        end
        if j<bih && i<biv
           tsum = tsum - Y2(c+oset2)/r2;
        end
        if j<bih && i>1 
           tsum = tsum - Y2(c+oset3)/r2;
        end

        G(c) = G(c) + beta*tsum;
    end
end


%----------------------------------------------------------------------
function [grad] = hubergradx(x, alpha)

    grad = max(min(x, alpha), -alpha);
    