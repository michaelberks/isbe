function C = mb_normxcorr2(template, image_in, weights, mask)
% mb version of the code below that allows template to have non-square
% shape, with missing values specified by NaNs
% In addition C is returned as 'valid' shape as opposed to 'full'
%NORMXCORR2 Normalized two-dimensional cross-correlation.
%   C = NORMXCORR2(TEMPLATE,A) computes the normalized cross-correlation of
%   matrices TEMPLATE and A. The matrix A must be larger than the matrix
%   TEMPLATE for the normalization to be meaningful. The values of TEMPLATE
%   cannot all be the same. The resulting matrix C contains correlation
%   coefficients and its values may range from -1.0 to 1.0.
%
%   Class Support
%   -------------
%   The input matrices can be numeric. The output matrix C is double.
%
%   Example
%   -------
%   template = .2*ones(11); % make light gray plus on dark gray background
%   template(6,3:9) = .6;   
%   template(3:9,6) = .6;
%   BW = template > 0.5;      % make white plus on black background
%   figure, imshow(BW), figure, imshow(template)
%   % make new image that offsets the template
%   offsetTemplate = .2*ones(21); 
%   offset = [3 5];  % shift by 3 rows, 5 columns
%   offsetTemplate( (1:size(template,1))+offset(1),...
%                   (1:size(template,2))+offset(2) ) = template;
%   figure, imshow(offsetTemplate)
%   
%   % cross-correlate BW and offsetTemplate to recover offset  
%   cc = normxcorr2(BW,offsetTemplate); 
%   [max_cc, imax] = max(abs(cc(:)));
%   [ypeak, xpeak] = ind2sub(size(cc),imax(1));
%   corr_offset = [ (ypeak-size(template,1)) (xpeak-size(template,2)) ];
%   isequal(corr_offset,offset) % 1 means offset was recovered
%
%  See also CORRCOEF.

%   Copyright 1993-2006 The MathWorks, Inc.
%   $Revision: 1.12.4.8 $  $Date: 2006/05/24 03:31:23 $

%   Input-output specs
%   ------------------
%   T:    2-D, real, full matrix
%         logical, uint8, uint16, or double
%         no NaNs, no Infs
%         prod(size(T)) >= 2
%         std(T(:))~=0
%
%   A:    2-D, real, full matrix
%         logical, uint8, uint16, or double
%         no NaNs, no Infs
%         size(A,1) >= size(T,1)
%         size(A,2) >= size(T,2)
%
%   C:    double
%

%   We normalize the cross correlation to get correlation coefficients using the
%   definition of Haralick and Shapiro, Volume II (p. 317), generalized to
%   two-dimensions. 
%
%   Lewis explicitly defines the normalized cross-correlation in two-dimensions
%   in this paper (equation 2):
%   
%      "Fast Normalized Cross-Correlation", by J. P. Lewis, Industrial Light & Magic.
%      http://www.idiom.com/~zilla/Papers/nvisionInterface/nip.html
%
%   Our technical reference document on NORMXCORR2 shows how to get from
%   equation 2 of the Lewis paper to the code below.

template = double(template);
nan_idx = isnan(template);


image_in = double(image_in);
if exist('mask', 'var') && ~isempty(mask)
    image_in(~mask) = 0;
end

%Convolve image with weights to compute local standard deviation of image
if ~exist('weights', 'var') || isempty(weights)
    weights = ones(size(template));
    weights(nan_idx) = 0;
end
weights = double(weights);
local_sum_A = do_correlation(image_in, weights);
local_sum_A2 = do_correlation(image_in.^2, weights);
    
%need to do this before we remove nan values


%Compute standard deviation of template. If we have a mask as input, this
%won't be constant over the image
if exist('mask', 'var') && ~isempty(mask)
    %Template size is not necessarily the same everywhere in the image
    mask = double(mask);
    template(nan_idx) = 0;
    size_T = do_correlation(double(mask), ones(size(template)));
    local_sum_T = do_correlation(mask, template);
    local_sum_T2 = do_correlation(mask, template.^2);
    diff_T_sums = local_sum_T2 - (local_sum_T.^2)./size_T;
    denom_T = sqrt(max(diff_T_sums,0) ); 

else
    size_T = sum(~nan_idx(:));
    denom_T = sqrt(size_T)*naNstd(template(:), 1);
    template(nan_idx) = 0;
    local_sum_T = sum(template(:));
end

%Convolve image with template
xcorr_TA = do_correlation(image_in, template);

% Note: diff_local_sums should be nonnegative, but may have negative
% values due to round off errors. Below, we use max to ensure the
% radicand is nonnegative.
diff_A_sums = local_sum_A2 - (local_sum_A.^2)./size_T;
denom_A = sqrt( max(diff_A_sums,0) ); 
denom = denom_T.*denom_A;

%numerator = (xcorr_TA/size_T) - (local_sum_A*naNsum(template(:))/size_T^2);
numerator = xcorr_TA - (local_sum_A.*local_sum_T)./size_T;

% We know denom_T~=0 from input parsing;
% so denom is only zero where denom_A is zero, and in 
% these locations, C is also zero.
tol = sqrt( eps( max(abs(denom(:)))) );
C = zeros(size(numerator));
i_nonzero = denom>tol;
C(i_nonzero) = numerator(i_nonzero) ./ denom(i_nonzero);
C(C>1) = 1;
C(C<-1)=-1;
if exist('mask', 'var') && ~isempty(mask)
    C(~mask) = 0;
end

%-------------------------------
% Function  do_correlation
%
function local_sum_A = do_correlation(A, shape)

% We thank Eli Horn for providing this code, used with his permission,
% to speed up the calculation of local sums. The algorithm depends on
% precomputing running sums as described in "Fast Normalized
% Cross-Correlation", by J. P. Lewis, Industrial Light & Magic.
% http://www.idiom.com/~zilla/Papers/nvisionInterface/nip.html

local_sum_A = conv2(A, rot90(shape, 2), 'same');

