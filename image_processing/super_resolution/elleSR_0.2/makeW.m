function [W] = makeW(ni_hi, nj_hi, ni_lo, nj_lo, theta, gamma)
%MAKEW computes the down-projected low-resolution image, from high-res image. Useful when using the generative model in the forward direction
%   [W, Y, theta] = makeW(v_hi, h_hi, lo_data)
%
%  *
%  * INPUTS:
%  *  ni_hi, nj_hi: vertical and horizontal sizes for the high-res image. 
%  *  ni_lo: a 1xK vector of lo-res image heights.
%  *  nj_lo: a 1xK vector of lo-res image widths.
%  *  theta: a 10xK matrix of geometric and photometric registration values.
%  *  gamma: a 1xK vector of noise std.
%  *
%  * OUTPUTS:
%  *  W: the super-res system matrix. Note that this is the TRANSPOSE of most of my papers' "W" matrices, since it's
%  *     rather easier to compute this way.
stdlim = 5; % controls how many standard deviations I truncate the PSF after. */

%Get the number of images and total number of lo-res pixels
K = length(ni_lo);
KM = sum(ni_lo .* nj_lo);

% One more loop through the inputs to start reading stuff into outputs. */
nzmax = 0;
for kj = 1:K %
    % this is pi*r^2*zoom*[number of pixels], where r=3*gamma. */
    nzmax = nzmax + ceil(pi*16*gamma(1, kj)*( theta(5,kj)*theta(1,kj) - theta(4,kj)*theta(2,kj) ))*ni_lo(kj)*nj_lo(kj); 
end
nzmax = max(nzmax, 5*10e6);

[W] = Wgauss(theta, gamma, ni_hi, nj_hi, nzmax, K, ni_lo, nj_lo, KM, stdlim);

%void Wgauss(mxArray **W, double *La, double *Lb, const double *theta, const int v_hi, const int h_hi, const int nzmax_in, const int K, const int *v_lo, const int* h_lo, const int KM) {
function [W] = Wgauss(theta, gamma, ni_hi, nj_hi, nzmax, K, ni_lo, nj_lo, KM, stdlim) 

h = zeros(ni_hi, nj_hi); % This is enough space to store a nonsparse row of W as I compute it. */
W_ridx = zeros(nzmax, 1);
W_cidx = zeros(nzmax, 1);
W_vals = zeros(nzmax, 1);

ypt = 0; %Tracks current point in the low-res image set
wpt = 1; %Tracks the current non-zero W entry

for kj = 1:K %(kj=0; kj<K; kj++) { % iterate over the K low-res images */
    
    thresh = stdlim*stdlim*gamma(kj);
    
    for j_lo = 1:nj_lo(kj) %(i_h_lo=0; i_h_lo < h_lo[kj]; i_h_lo++) { % iterate over low-res image's rows. */
        for i_lo = 1:ni_lo(kj) %(i_v_lo=0; i_v_lo < v_lo[kj]; i_v_lo++) { % scan down each column in turn. */
 
            %Column index for this pixel in Y(k)
            W_ci = ypt + (j_lo-1)*ni_lo(kj) + i_lo;
            
            % This call maps the location to the HR frame and finds the params of the affine approximation to the PSF under projection. */
            [i_c, j_c, b11, b12, b22, delta_i, delta_j] =  getGaussParams(theta(:,kj), gamma(kj), i_lo, j_lo, stdlim);
            
            if (j_c >=0 && j_c<nj_hi && i_c>=0 && i_c<ni_hi) % this means the mapped low-res pixel lands somewhere in the superimage. */
                j0 = max(floor(j_c-delta_j), 1);
                i0 = max(floor(i_c-delta_i), 1);
                j1 = min(ceil(j_c+delta_j), nj_hi);
                i1 = min(ceil(i_c+delta_i), ni_hi); % Find the right box in the HR image to scan over. */
                norm_val = 0;
                
                % Find all the values at pixel locations (and for a 1-pixel border around region). */
                for j_hi = j0:j1
                    for i_hi = i0:i1
                        di = i_hi-i_c;
                        dj = j_hi-j_c;
                        tval = b11*dj*dj + b12*dj*di + b22*di*di;
                        if (tval<thresh)
                            h(i_hi, j_hi) = exp(-0.5*tval); % This evaluates the Gaussian kernel value. */
                            norm_val = norm_val + h(i_hi, j_hi); % I'll need this value later for normalization. */
                        else
                            h(i_hi, j_hi) = 0; 
                        end
                    end
                end
                if (norm_val > 0.0000001) % if all the values are low, I'll flag the pixel as unsupported rather than assume good values. */
                    for j_hi = j0:j1 % iterate back across the patch, storing the normalized values in W. */
                        for i_hi = i0:i1
                            if h(i_hi, j_hi)>0
                                %Column index for this pixel in hi-es image X
                                W_cidx(wpt) = (j_hi-1)*ni_hi+i_hi;
                                
                                %Store row as computed above
                                W_ridx(wpt) = W_ci;
                                
                                % data element for W matrix (normalized). */
                                W_vals(wpt) = h(i_hi,j_hi) / norm_val;
                                
                                %Increment wpt
                                wpt = wpt + 1;
                            end
                        end
                    end
                else
                    sprintf('WARNING! Insufficient support for psf!\n'); % This happens if "myscale" is too low. */
                end
            end
        end
    end
    ypt = ypt + nj_lo(kj)*ni_lo(kj);
end

% %Remove any indices and values we didn't use
% W_ridx(wpt:end) = [];
% W_cidx(wpt:end) = [];
% W_vals(wpt:end) = [];

%Build the sparse system matrix W from the indices and values
W = sparse(W_ridx(1:wpt-1), W_cidx(1:wpt-1), W_vals(1:wpt-1), KM, ni_hi*nj_hi);

% All done! */

%-------------------------------------------------------------------------------------
%getLambdaGauss(const double *M_k, const double i_v_lo, const double i_h_lo, double *nuv, 
%   double *nuh, double *b11, double *b12, double *b22, double *deltav, double *deltah)
function [i_c, j_c, b11, b12, b22, delta_i, delta_j] =  getGaussParams(theta_k, gamma_k, i_lo, j_lo, stdlim)

denom = j_lo*theta_k(7) + i_lo*theta_k(8) + 1;
j_c = (j_lo*theta_k(1) + i_lo*theta_k(2)  + theta_k(3))/denom; %
i_c = (j_lo*theta_k(4) + i_lo*theta_k(5)  + theta_k(6))/denom; %
denom = denom * denom;

% H = [h11,h12;h21,h22] is the Hessian of the transform given by M_k, evaluated at
%   the point (i_v_lo,i_h_lo) maps to under that transform. */
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
delta_j = 2*stdlim*gam*sqrt(b22)/denom; % vertical extent of the new psf kernel */
return;
