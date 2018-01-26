function [avg_im, x_mask] = getAvim(ni_hi, nj_hi, ni_lo, nj_lo, Y, theta, gamma, Y_mask)
%  [avim,msk,theta] = getAvim(ni_hi,nj_hi,o);
%  *
%  * INPUTS:
%  *  ni_hi, nj_hi: vertical and horizontal sizes for the average image. 
%  *  o: the super-res problem datastructure, which must be a MATLAB struct with the following entries:
%  *  theta: a 10xK matrix of geometric registration values.
%  *    o(i).H is a 3x3 double representing the homography taking points in the i^th LR image into the SR image.
%  *    o(i).la is the multiplicative photometric parameter for image i.
%  *    o(i).lb is the additive photometric parameter for image i.
%  *    o(i).g is the sigma for the Gaussian PSF associated with image i.
%  *    o(i).im is the i^th image. W doesn't depend on image pixels, but the image size is important.
%  *
%  * OUTPUTS:
%  *  avim: the average image (which will be of size ni_hi-by-nj_hi).
%  *  msk: a double matrix of size ni_hi-by-nj_hi indicating which pixels in avim were estimated.
%  *    -- average image pixels with negligible support from any of the PSFs of the low-res
%  *       image pixels are returned as zeros.

if ~exist('Y_mask', 'var') || isempty(Y_mask)
    Y_mask = true(size(Y));
end

stdlim = 5; % controls how many standard deviations I truncate the PSF after. */
% Now deal with input prhs[2], which is the K-element structure, o. */

K = size(theta, 2);

if length(ni_lo) == 1
    ni_lo = repmat(ni_lo, K, 1);
end
if length(nj_lo) == 1
    nj_lo = repmat(nj_lo, K, 1);
end

avg_im = zeros(ni_hi, nj_hi);
x_mask = zeros(ni_hi, nj_hi);
h = zeros(ni_hi, nj_hi);

curr_pt = 0;
for kj = 1:K %(kj=0; kj<K; kj++) { % iterate over the K low-res images */
    
    %Reshape the lo-res image from the vectorised Y
    n_pts = ni_lo(kj)*nj_lo(kj);
    idx = curr_pt + (1:n_pts);
    curr_pt = curr_pt + n_pts;
    im = reshape(Y(idx), ni_lo(kj), nj_lo(kj));
    im_mask = reshape(Y_mask(idx), ni_lo(kj), nj_lo(kj));
    
    for j_lo = 1:nj_lo(kj) % iterate over low-res image's rows. */
        for i_lo = 1:ni_lo(kj) % scan down each column in turn. */
            
            if ~im_mask(i_lo, j_lo)
                continue;
            end

            % This call maps the location to the HR frame and find the params of the affine approximation to the PSF under projection. */
            [i_c, j_c, b11, b12, b22, delta_i, delta_h] =  getGaussParams(theta(:,kj), gamma(:,kj), i_lo, j_lo, stdlim);

            if (j_c >=1 && j_c<nj_hi && i_c>=0 && i_c<ni_hi) % this means the mapped low-res pixel lands somewhere in the superimage. */
                j0 = max(floor(j_c-delta_h), 1);
                i0 = max(floor(i_c-delta_i), 1);
                j1 = min(ceil(j_c+delta_h), nj_hi);
                i1 = min(ceil(i_c+delta_i), ni_hi); % Find the right box in the HR image to scan over. */
                norm_val = 0;
                % Find all the values at pixel locations (and for a 1-pixel border around region). */
                for j_hi = j0:j1 %(ibih = j0; ibih < j1; ibih++) {
                    for i_hi = i0:i1 %(ibiv = i0; ibiv < i1; ibiv++) {
                        nv = i_hi-i_c;
                        nh = j_hi-j_c;
                        h(i_hi, j_hi) = exp(-0.5*(b11*nh*nh + b12*nh*nv + b22*nv*nv));
                        norm_val = norm_val + h(i_hi, j_hi);
                    end
                end
                if (norm_val>0.0000001) 
                    for j_hi = j0:j1 %(ibih = j0; ibih < j1; ibih++) {
                        for i_hi = i0:i1 %(ibiv = i0; ibiv < i1; ibiv++) {
                            h(i_hi, j_hi) = h(i_hi, j_hi)/norm_val;
                            avg_im(i_hi, j_hi) = avg_im(i_hi, j_hi) + im(i_lo, j_lo)*h(i_hi, j_hi);
                            x_mask(i_hi, j_hi) = x_mask(i_hi, j_hi) + h(i_hi, j_hi);
                        end
                    end
                end
            end
        end
    end
end

for j_hi = 1:nj_hi
    for i_hi = 1:ni_hi
        if (x_mask(i_hi, j_hi)>0.0000001)
            avg_im(i_hi, j_hi) = avg_im(i_hi, j_hi)/x_mask(i_hi, j_hi);
            x_mask(i_hi, j_hi) = 0;
        else 
            avg_im(i_hi, j_hi) = 0;
            x_mask(i_hi, j_hi) = 1;
        end
    end
end

%-------------------------------------------------------------------------------------
%getLambdaGauss(const double *theta_k, const double ismv, const double ismh, double *i_c, 
%   double *j_c, double *b11, double *b12, double *b22, double *delta_i, double *delta_h)
function [i_c, j_c, b11, b12, b22, delta_i, delta_h] =  getGaussParams(theta_k, gamma_k, i_lo, j_lo, stdlim)

denom = j_lo*theta_k(7) + i_lo*theta_k(8) + 1;
j_c = (j_lo*theta_k(1) + i_lo*theta_k(2)  + theta_k(3))/denom; %
i_c = (j_lo*theta_k(4) + i_lo*theta_k(5)  + theta_k(6))/denom; %
denom = denom * denom;

% H = [h11,h12;h21,h22] is the Hessian of the transform given by theta_k, evaluated at
%   the point (ismv,ismh) maps to under that transform. */
h11 = ((theta_k(1)*theta_k(8) - theta_k(2)*theta_k(7))*i_lo + theta_k(1) - theta_k(3)*theta_k(7))/denom;
h12 = ((theta_k(2)*theta_k(7) - theta_k(1)*theta_k(8))*j_lo + theta_k(2) - theta_k(3)*theta_k(8))/denom;
h21 = ((theta_k(4)*theta_k(8) - theta_k(5)*theta_k(7))*i_lo + theta_k(4) - theta_k(6)*theta_k(7))/denom;
h22 = ((theta_k(5)*theta_k(7) - theta_k(4)*theta_k(8))*j_lo + theta_k(5) - theta_k(6)*theta_k(8))/denom;

% I also know that the covariance of the psf will be H*Sig*H', where sig was the
%   original covariance in the low-res image, which has a variance of sqgam. */
detH = h11*h22 - h12*h21;
detH = 1/(gamma_k*(detH*detH));
b11 = detH*(h21*h21+h22*h22); % inv(H*Sig*H')(1,1) */
b12 = -2*detH*(h11*h21+h12*h22);% inv(H*Sig*H')(1,2)+inv(H*Sig*H')(2,1) */
b22 = detH*(h11*h11+h12*h12); % inv(H*Sig*H')(2,2) */

% Now find the greatest h and v extend of this PSF that's within 4*gam in the original (low-res) PSF.
% * Derive this by expressing x in terms of y, using the quadratic formula, and finding the value of y^2 necessary to set the sqrt(b^2-4ac) part to zero. */
denom = sqrt(4*b11*b22 - b12*b12);
gam = sqrt(gamma_k);
delta_i = 2*stdlim*gam*sqrt(b11)/denom; % vertical extent of the new psf kernel */
delta_h = 2*stdlim*gam*sqrt(b22)/denom; % vertical extent of the new psf kernel */
