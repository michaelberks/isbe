function [variance] = imf_sigma(input_im, mask, roi)
%IMF_SIGMA *Insert a one line summary here*
%   [variance] = imf_sigma(input_im, mask, roi)
%
% Inputs:
%      input_im - *Insert description of input variable here*
%
%      mask - *Insert description of input variable here*
%
%      roi - *Insert description of input variable here*
%
%
% Outputs:
%      variance - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 29-Jun-2017
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

if ~exist('mask', 'var') || isempty(mask)
    mask = true(size(input_im));
end
if ~exist('roi', 'var') || isempty(roi)
    roi = [1 size(input_im,1) 1 size(input_im,2)];
end
mask([1:roi(1) roi(2)+1:end],:) = 0;
mask(:,[1:roi(3) roi(4)+1:end]) = 0;
    
hist_struc = initialise_histogram(input_im(mask), [], [], 50);
[max_val, max_x] = max_hist(hist_struc);

%Find the first bin either side of the peak with count less than half the
%peak
fwhm = zeros(2,1);
for i = 1:2
    num = 1;
    height = max_val;
    while (height > 0.5*max_val)
        height = hist_val(hist_struc, max_x+power(-1.0,i)*num*hist_struc.xincr);
        num = num+1;
    end
    fwhm(i) = max_x + power(-1.0,i)*num*(hist_struc.xincr);
end

%Create a new range that...
range = 3.0*(fwhm(2)-fwhm(1))/2.35;

% Recompute fitted histograms with lower range?
hist_struc = initialise_histogram(input_im(mask), [], [-range range], 50);

%Get hist value at positive and negative tails
maxtailp = hist_struc.xcounts(end);
maxtailm = hist_struc.xcounts(1);

%Compute background as average of positive and negative tails
bkg =  (maxtailp+maxtailm)/2;  
hist_struc = add_histogram_data(hist_struc, hist_struc.xbins, -bkg*ones(size(hist_struc.xbins)));

change = 0;
for step = hist_struc.xbins
    change = change + hist_val(hist_struc,step)*abs(step-hist_struc.mean);
end

% /* include correction for truncation at 1.5 sigma */
% /*
% variance = 2.0*fabs(hist_struc.mean2-hist_struc.mean*hist_struc.mean);
% */
if (hist_struc.contents > eps) 
    variance = 1.035*change/hist_struc.contents;
else
    variance = 0.0;
end

%***None of the below seem to do anything but were in the Tina code***
% %Compute cut as ratio of max peak to background dist
% %Get new max val
% [max_val, max_x] = max_hist(hist_struc);
% 
% %Compute tail cuttoff?
% htail = 0.0333/(range*sqrt(2*3.14159));
% cut = (maxval-bkg)*0.325;
% if (abs(hist_struc.mean2 - hist_struc.mean*hist_struc.mean) > eps)
%      graph_hfit(imcalc_graph_tv_get(), hist);
% end
