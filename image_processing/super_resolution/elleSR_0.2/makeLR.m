function [im,noise,orig] = makeLR(X, ni_lo, nj_lo, theta, lambda, sigma, gamma, mask)
%MAKELR computes the down-projected low-resolution image, from high-res image. Useful when using the generative model in the forward direction
%   [] = makeLR()
%
%  * INPUTS:
%  *  X: high-res image. 
%  *  o: the super-res problem datastructure, which must be a MATLAB struct with the following entries:
%  *    o.H is a 3x3 double representing the homography taking points in the LR image into the SR image.
%  *    o.la is the multiplicative photometric parameter.
%  *    o.lb is the additive photometric parameter.
%  *    o.g is the sigma for the Gaussian PSF.
%  *    o.n is the sigma for the additive Gaussian iid noise on each pixel.
%  *    o.v vertical size for the image.
%  *    o.h horizontal size for the image.
%  *  NOTE: if o is larger than 1x1, only the first part is used (i.e. this function doens't 
%  *        currently produce a stack of images, just one at a time).
%  *
%  * OUTPUTS:
%  *  im: the generated low-resolution image, including image noise.
%  *  noise: the noise image.
%  *  orig: the low-resolution image without the added noise.
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

if ~exist('mask', 'var') || isempty(mask)
    mask = true(ni_lo, nj_lo);
elseif size(mask,1) ~= ni_lo
    mask = reshape(mask, ni_lo, nj_lo);
end 

% Grab inputs */
[ni_hi, nj_hi] =  size(X);

orig = getYgauss(ni_hi, nj_hi, ni_lo, nj_lo, gamma, theta, lambda, X, mask);
noise = randn(ni_lo, nj_lo)*sigma;
im = orig + noise;

%--------------------------------------------------------------------------
%void getLambdaGauss(const double *M_k, const double sqgam, const double ismv, 
%   const double ismh, double *lambda, double *nuv, double *nuh, double *b11, double *b12, double *b22)
function [phi, i_c, j_c, b11, b12, b22] = getGaussParams(theta, sqgam, i_lo,  j_lo)

kdenom = j_lo*theta(7) + i_lo*theta(8) + 1;
j_c = (j_lo*theta(1) + i_lo*theta(2)  + theta(3))/kdenom; %
i_c = (j_lo*theta(4) + i_lo*theta(5)  + theta(6))/kdenom; %
kdenom = kdenom * kdenom;

% H = [h11,h12;h21,h22] is the Hessian of the transform given by M_k, evaluated at
%   the point (ismv,ismh) maps to under that transform.
h11 = ((theta(1)*theta(8) - theta(2)*theta(7))*i_lo + theta(1) - theta(3)*theta(7))/kdenom;
h12 = ((theta(2)*theta(7) - theta(1)*theta(8))*j_lo + theta(2) - theta(3)*theta(8))/kdenom;
h21 = ((theta(4)*theta(8) - theta(5)*theta(7))*i_lo + theta(4) - theta(6)*theta(7))/kdenom;
h22 = ((theta(5)*theta(7) - theta(4)*theta(8))*j_lo + theta(5) - theta(6)*theta(8))/kdenom;

% I also know that the covariance of the psf will be H*Sig*H', where sig was the
%   original covariance in the low-res image, which has a variance of sqgam. */
detH = h11*h22 - h12*h21;
detH = 1/(sqgam*(detH*detH));
b11 = detH*(h21*h21+h22*h22); % inv(H*Sig*H')(1,1) */
b12 = -2*detH*(h11*h21+h12*h22); % inv(H*Sig*H')(1,2)+inv(H*Sig*H')(2,1) */
b22 = detH*(h11*h11+h12*h12); % inv(H*Sig*H')(2,2) */

% The scaling along the semimajor axis of the new covariance ellipse will be given by the sqrt of the
%   bigger eigenvalue of inv(Sig). I'm going to call that scaling 'phi' in here. */
phi = b11 + b22;
phi = 0.5*(phi - sqrt(phi*phi - 4*(b11*b22 - 0.25*b12*b12)));
phi = 3/sqrt(phi) + 1; % Multiplying by 3 means everything with P>tiny lies within this distance of the transformed mu. */

%--------------------------------------------------------------------------
% getY gets the down-projected superim Yd. Useful when using the generative model in the forward direction. */
% void getYgauss(double *Yd, const int v_hi, const int h_hi, const int v_lo, const int h_lo, const double sqgam, const double* Hph, const double* x) {

function [Yd] = getYgauss(ni_hi, nj_hi, ni_lo, nj_lo, gamma_2, theta, lambda, x, mask)

Yd = zeros(ni_lo, nj_lo);

for j_lo = 1:nj_lo %
    for i_lo = 1:ni_lo %
        
        if ~mask(i_lo, j_lo)
            continue;
        end
            
        [phi, i_c, j_c, b11, b12, b22] = getGaussParams(theta, gamma_2, i_lo, j_lo);

        if ( j_c>0 && j_c<nj_hi && i_c>0 && i_c<ni_hi) %% this means the mapped low-res pixel lands somewhere in the superimage. */
            j0 = max(floor(j_c-phi), 1);
            i0 = max(floor(i_c-phi), 1);
            j1 = min(ceil(j_c+phi), nj_hi);
            i1 = min(ceil(i_c+phi), ni_hi);
            norm_val = 0;

            % Find all the values at pixel locations (and for a 1-pixel border around region). */
            for j_hi = j0:j1
                for i_hi = i0:i1
                    di = i_hi-i_c;
                    dj = j_hi-j_c;
                    w_val = exp(-0.5*(b11*dj*dj + b12*dj*di + b22*di*di));
                    norm_val = norm_val + w_val;
                    Yd(i_lo, j_lo) = Yd(i_lo, j_lo) + w_val*x(i_hi, j_hi);
                end
            end
            if (norm_val > 0.0000001) 
                Yd(i_lo, j_lo) = (lambda(1)*(Yd(i_lo, j_lo)/norm_val) + lambda(2)); 

            else
                display('WARNING! Insufficient support for psf!\n'); 
                sprintf('( %i , %i ) -> ( %g , %g ); lambda = %g \n',i_lo,j_lo,i_c,j_c,phi);
                sprintf('v %i to %i; h %i to %i \n',i0,i1,j0,j1);
                sprintf('b = ( %g , %g , %g ] \n',b11,b12,b22);
                return;
            end
        end
    end
end
