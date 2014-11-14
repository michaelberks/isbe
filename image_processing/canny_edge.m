function [final_edges ori ax ay] = canny_edge(image_in, thresh, sigma, percent_not_edges, thresh_ratio, ignore_map, if_plot)

    if nargin < 2
        %Check if we've been given thresholds
        thresh = [];
    end
    if nargin < 3
        %Check if sigma has been set
        sigma = 1;
    end
    if nargin < 4
        percent_not_edges = .99;
    end
    if nargin < 5
        thresh_ratio = .7;
    end
    
    image_in = double(image_in); 
    [m n] = size(image_in);

        
    if nargin < 6
        ignore_map = false(m,n);
    end
    if nargin < 7
        if_plot = false;
    end
    
    % Magic numbers
    GaussianDieOff = .0001;  

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
    aSmooth=imfilter(image_in,gau,'conv','replicate');   % run the filter across rows
    aSmooth=imfilter(aSmooth,gau','conv','replicate'); % and then across columns

    %apply directional derivatives
    ax = imfilter(aSmooth, dgau2D, 'conv','replicate');
    ay = imfilter(aSmooth, dgau2D', 'conv','replicate');
    
    %Get the magnitude of the filter response
    mag = sqrt((ax.*ax) + (ay.*ay));
    magmax = max(mag(:));
    if magmax > 0
        mag = mag / magmax;   % normalize
    end
    
    %Now do the maximal suppression
    %compute orientations at each point
    ori = atan(ay./ax);
    %We may have some NaNs where ax and ay are both zero, so set these to
    %zero
    ori(isnan(ori)) = 0;
    ori_c = cos(ori);
    ori_s = sin(ori);

    %compute (x, y) coordinates along +/- unit normal vectors
    xx = repmat(1:n, m, 1);
    yy = repmat((1:m)', 1, n);
    
    x_interp1 = xx + ori_c;
    x_interp2 = xx - ori_c;

    y_interp1 = yy + ori_s;
    y_interp2 = yy - ori_s;

    %linearly interpolate values at normal coordinates
    z1 = interp2(mag, x_interp1, y_interp1);
    z2 = interp2(mag, x_interp2, y_interp2);

    %discard any point from im that is smaller than either of its two normals 
    discard = (mag <= z1) | (mag <= z2);
    nms_mag = mag;
    nms_mag(discard) = 0;
    
    %Now apply the ignore map to throw away edges from a region we don't
    %want
    nms_mag(ignore_map) = 0;
    nms_mag(ignore_map) = 0;

    % Select the thresholds
    if isempty(thresh) 
        counts = imhist(nms_mag, 64);
        high_thresh = ...
            find(cumsum(counts) > percent_not_edges*m*n, 1,'first') / 64;
        low_thresh = thresh_ratio*high_thresh;

    elseif length(thresh)==1
        high_thresh = thresh;
        if thresh >= 1
            eid = sprintf('Images:%s:thresholdMustBeLessThanOne', mfilename);
            msg = 'The threshold must be less than 1.'; 
            error(eid,'%s',msg);
        end
        low_thresh = thresh_ratio*thresh;

    elseif length(thresh) == 2
        low_thresh = thresh(1);
        high_thresh = thresh(2);
        
        if (low_thresh >= high_thresh) || (high_thresh >= 1)
            eid = sprintf('Images:%s:thresholdOutOfRange', mfilename);
            msg = 'Thresh must be [low high], where low < high < 1.'; 
            error(eid,'%s',msg);
        end
    end
    
    %Now do the edge hysterisis
    %Compute the strong edges and the weak edges
    weak_edges = nms_mag > low_thresh;
    strong_edges = nms_mag > high_thresh;
    
    %Now find weak edges that are 8-connected to a strong edge
    if any(strong_edges(:))
        [rstrong cstrong] = find(strong_edges);
        combined_edges = bwselect(weak_edges, cstrong, rstrong, 8);
    else
        combined_edges = weak_edges;
    end
    
    %finally, thin the combined edges
    final_edges = bwmorph(combined_edges, 'thin', 1);
    
    if if_plot
        figure; 
        subplot(2,3,1); imagesc(mag); axis image; colormap(gray(256)); title('Original magnitudes');
        subplot(2,3,2); imagesc(nms_mag); axis image; colormap(gray(256)); title('Suppressed magnitudes');
        subplot(2,3,3); imagesc(weak_edges); axis image; colormap(gray(256)); title('Weak edges');
        subplot(2,3,4); imagesc(strong_edges); axis image; colormap(gray(256)); title('Strong edges');
        subplot(2,3,5); imagesc(combined_edges); axis image; colormap(gray(256)); title('Combined edges');
        subplot(2,3,6); imagesc(final_edges); axis image; colormap(gray(256)); title('Final edges');
    end