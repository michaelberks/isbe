function e = canny_edge_horizontal(image_in, thresh, sigma)

    if nargin < 2
        thresh = [];
    end
    if nargin < 3
        sigma = 1;
    end
    
    [m n] = size(image_in);
    
    % Magic numbers
    GaussianDieOff = .0001;  
    PercentOfPixelsNotEdges = .7; % Used for selecting thresholds
    ThresholdRatio = .4;          % Low thresh is this fraction of the high.

    % Design the filters - a gaussian and its derivative

    pw = 1:30; % possible widths
    ssq = sigma^2;
    width = find(exp(-(pw.*pw)/(2*ssq))>GaussianDieOff,1,'last');
    if isempty(width)
        width = 1;  % the user entered a really small sigma
    end

    t = (-width:width);
    gau = exp(-(t.*t)/(2*ssq))/(2*pi*ssq);     % the gaussian 1D filter

    % Find the directional derivative of 2D Gaussian (along X-axis)
    % Since the result is symmetric along X, we can get the derivative along
    % Y-axis simply by transposing the result for X direction.
    [x,y]=meshgrid(-width:width,-width:width);
    dgau2D=-x.*exp(-(x.*x+y.*y)/(2*ssq))/(pi*ssq);

    % Convolve the filters with the image in each direction
    % The canny edge detector first requires convolution with
    % 2D gaussian, and then with the derivitave of a gaussian.
    % Since gaussian filter is separable, for smoothing, we can use 
    % two 1D convolutions in order to achieve the effect of convolving
    % with 2D Gaussian.  We convolve along rows and then columns.

    %smooth the image out
    aSmooth=imfilter(image_in,gau,'conv','replicate');   % run the filter accross rows
    aSmooth=imfilter(aSmooth,gau','conv','replicate'); % and then accross columns

    %apply directional derivatives
    ax = imfilter(aSmooth, dgau2D, 'conv','replicate');
    ay = imfilter(aSmooth, dgau2D', 'conv','replicate');

    %mag = sqrt((ax.*ax) + (ay.*ay));
    mag = abs(ax);
    magmax = max(mag(:));
    if magmax > 0
        mag = mag / magmax;   % normalize
    end

    % Select the thresholds
    if isempty(thresh) 
        counts=imhist(mag, 64);
        highThresh = find(cumsum(counts) > PercentOfPixelsNotEdges*m*n,...
        1,'first') / 64;
        lowThresh = ThresholdRatio*highThresh;
        thresh = [lowThresh highThresh];
    elseif length(thresh)==1
        highThresh = thresh;
        if thresh>=1
            eid = sprintf('Images:%s:thresholdMustBeLessThanOne', mfilename);
            msg = 'The threshold must be less than 1.'; 
            error(eid,'%s',msg);
        end
        lowThresh = ThresholdRatio*thresh;
        thresh = [lowThresh highThresh];
    elseif length(thresh)==2
        lowThresh = thresh(1);
        highThresh = thresh(2);
        if (lowThresh >= highThresh) || (highThresh >= 1)
            eid = sprintf('Images:%s:thresholdOutOfRange', mfilename);
            msg = 'Thresh must be [low high], where low < high < 1.'; 
            error(eid,'%s',msg);
        end
    end

    % The next step is to do the non-maximum supression.  
    % We will accrue indices which specify ON pixels in strong edgemap
    % The array e will become the weak edge map.
    idxStrong = [];
    e = zeros(size(mag));
    for dir = [1 4]
        idxLocalMax = cannyFindLocalMaxima(dir,ax,ay,mag);
        idxWeak = idxLocalMax(mag(idxLocalMax) > lowThresh);
        e(idxWeak)=1;
        idxStrong = [idxStrong; idxWeak(mag(idxWeak) > highThresh)]; %#ok
    end

    if ~isempty(idxStrong) % result is all zeros if idxStrong is empty
        rstrong = rem(idxStrong-1, m)+1;
        cstrong = floor((idxStrong-1)/m)+1;
        e = bwselect(e, cstrong, rstrong, 8);
        e = bwmorph(e, 'thin', 1);  % Thin double (or triple) pixel wide contours
    end

function idxLocalMax = cannyFindLocalMaxima(direction,ix,iy,mag)
%
% This sub-function helps with the non-maximum supression in the Canny
% edge detector.  The input parameters are:
% 
%   direction - the index of which direction the gradient is pointing, 
%               read from the diagram below. direction is 1, 2, 3, or 4.
%   ix        - input image filtered by derivative of gaussian along x 
%   iy        - input image filtered by derivative of gaussian along y
%   mag       - the gradient magnitude image
%
%    there are 4 cases:
%
%                         The X marks the pixel in question, and each
%         3     2         of the quadrants for the gradient vector
%       O----0----0       fall into two cases, divided by the 45 
%     4 |         | 1     degree line.  In one case the gradient
%       |         |       vector is more horizontal, and in the other
%       O    X    O       it is more vertical.  There are eight 
%       |         |       divisions, but for the non-maximum supression  
%    (1)|         |(4)    we are only worried about 4 of them since we 
%       O----O----O       use symmetric points about the center pixel.
%        (2)   (3)        


[m,n] = size(mag);

% Find the indices of all points whose gradient (specified by the 
% vector (ix,iy)) is going in the direction we're looking at.  

switch direction
 case 1
  idx = find((iy<=0 & ix>-iy)  | (iy>=0 & ix<-iy));
 case 2
  idx = find((ix>0 & -iy>=ix)  | (ix<0 & -iy<=ix));
 case 3
  idx = find((ix<=0 & ix>iy) | (ix>=0 & ix<iy));
 case 4
  idx = find((iy<0 & ix<=iy) | (iy>0 & ix>=iy));
end

% Exclude the exterior pixels
if ~isempty(idx)
  v = mod(idx,m);
  extIdx = find(v==1 | v==0 | idx<=m | (idx>(n-1)*m));
  idx(extIdx) = [];
end

ixv = ix(idx);  
iyv = iy(idx);   
gradmag = mag(idx);

% Do the linear interpolations for the interior pixels
switch direction
 case 1
  d = abs(iyv./ixv);
  gradmag1 = mag(idx+m).*(1-d) + mag(idx+m-1).*d; 
  gradmag2 = mag(idx-m).*(1-d) + mag(idx-m+1).*d; 
 case 2
  d = abs(ixv./iyv);
  gradmag1 = mag(idx-1).*(1-d) + mag(idx+m-1).*d; 
  gradmag2 = mag(idx+1).*(1-d) + mag(idx-m+1).*d; 
 case 3
  d = abs(ixv./iyv);
  gradmag1 = mag(idx-1).*(1-d) + mag(idx-m-1).*d; 
  gradmag2 = mag(idx+1).*(1-d) + mag(idx+m+1).*d; 
 case 4
  d = abs(iyv./ixv);
  gradmag1 = mag(idx-m).*(1-d) + mag(idx-m-1).*d; 
  gradmag2 = mag(idx+m).*(1-d) + mag(idx+m+1).*d; 
end
idxLocalMax = idx(gradmag>=gradmag1 & gradmag>=gradmag2); 