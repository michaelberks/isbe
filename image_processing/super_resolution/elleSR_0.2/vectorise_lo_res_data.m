function [Y, theta, lambda, gamma, ni_lo, nj_lo] = vectorise_lo_res_data(lo_data)
%MAKEW computes the down-projected low-resolution image, from high-res image. Useful when using the generative model in the forward direction
%   [W, Y, theta] = makeW(v_hi, h_hi, lo_data)
%
%  *
%  * INPUTS:
%  *  lo_data: the super-res problem datastructure, which must be a MATLAB struct with the following entries:
%  *    lo_data(i).H is a 3x3 double representing the homography taking points in the i^th LR image into the SR image.
%  *    lo_data(i).la is the multiplicative photometric parameter for image i.
%  *    lo_data(i).lb is the additive photometric parameter for image i.
%  *    lo_data(i).g is the sigma for the Gaussian PSF associated with image i.
%  *    lo_data(i).im is the i^th image. W doesn't depend on image pixels, but the image size is important.
%  *
%  * OUTPUTS:
%  *  W: the super-res system matrix. Note that this is the TRANSPOSE of most of my papers' "W" matrices, since it's
%  *     rather easier to compute this way.
%  *  Y: A KMx1 vector containing all the low-res pixel values, where "KM" is the total number of low-res pixels in the data.
%  *  theta: a 10xK matrix of geometric registration values.
%  *  gamma: a 1xK vector of geometric registration values.
%  *  v_lo: a 1xK vector of lo-res image heights.
%  *  h_lo: a 1xK vector of lo-res image widths.

% Check for proper number of arguments. */

% Now deal with input prhs[2], which is the K-element structure, lo_data. */
if ~isstruct(lo_data)
    error('Input must be a structure.');
end

K = length(lo_data);
ni_lo = zeros(K,1);
nj_lo = zeros(K,1);

% Check all the contents of the struct I need are actually there */
required_fieldnames = {'H', 'l1', 'l2', 'g', 'im'};
supplied_fieldnames = fieldnames(lo_data);
if length(intersect(required_fieldnames, supplied_fieldnames)) ~= length(required_fieldnames)
    for i_f = 1:length(required_fieldnames)
        if ~ismember(required_fieldnames{i_f}, supplied_fieldnames)
            sprintf('Input struct lo_data is missing field %s\n', required_fieldnames{i_f});
        end
    end
    error('Not all required fieldnames are present in lo_data');
end

%Get size of each low-res image and compute total number of low res pixels
for kj = 1:K %(kj=0;kj<K;kj++) {
    [ni_lo(kj), nj_lo(kj)] = size(lo_data(kj).im);
end
KM = sum(ni_lo.*nj_lo);

% Now create the output matrices. */
Y = zeros(KM,1); % Low-res pixels vectorised
theta = zeros(8,K);
lambda = zeros(2,K);
gamma = zeros(1, K);

% One more loop through the inputs to start reading stuff into outputs. */
ypt = 0;
for kj = 1:K %(kj=0;kj<K;kj++) {
    
    num_pixels = ni_lo(kj)*nj_lo(kj);
    idx = ypt + (1:num_pixels);
    
    Y(idx) = lo_data(kj).im(:);
    ypt = ypt + num_pixels;
    
    H = lo_data(kj).H;
    theta(1, kj) = H(1)/H(9); % I want theta to represent row-major rather than         */
    theta(2, kj) = H(4)/H(9); % column-major form, because that's what I            */
    theta(3, kj) = H(7)/H(9); % wrote the getLambdaGauss function to work with..    */
    theta(4, kj) = H(2)/H(9); %  .                                                  */
    theta(5, kj) = H(5)/H(9); %  .                                                  */
    theta(6, kj) = H(8)/H(9); %  .                                                  */
    theta(7, kj) = H(3)/H(9); %  .                                                  */
    theta(8, kj) = H(6)/H(9); % .. so these lines are doing a transpose.            */
    lambda(1, kj) = lo_data(kj).l1;
    lambda(2, kj) = lo_data(kj).l2;
    gamma(1, kj) = lo_data(kj).g^2;
end