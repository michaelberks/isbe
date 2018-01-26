function Dx = diff_x(X, ni, nj)
%DIFF_X
%   [] = diff_x(X, ni, nj)
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
Dx = zeros(length(X),1);
   
% huber_prior_gradx(G, X, ni, nj, nu) */
% adds the gradient of p(x) wrt x to the current values of the elements of G. */
r2 = sqrt(2);

for i = 1:ni %(i=0; i<ni; i++) {
    for j = 1:nj %(j=0; j<nj; j++) {
        c = (j-1)*ni+i;
        if i < ni 
            Dx(c,1) = -X(c)+X(c+1); %Horizontal grad
            if j<nj 
                Dx(c, 2) = (-X(c)+X(c+ni+1))/r2; %Vertical grad 
            end
        end
        if j<nj
            Dx(c, 3) = -X(c)+X(c+ni); %Diag L
            if i >= 2  
                Dx(c, 4) = (-X(c)+X(c+ni-1))/r2; %Diag R
            end
        end
    end
end