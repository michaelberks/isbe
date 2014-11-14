function [response orientation] = gabor_filter(varargin)

%
% GABOR_FILTER Apply Gabor filter to an image
%
% This function applies a bank of Gabor Filters to an image. The arguments
% tau and len are defined in Rangayyan and Ayres 2006
%
% GABOR_FILTER uses the U_PACKARGS interface function
% and so arguments to the function may be passed as name-value pairs
% or as a struct with fields with corresponding names. The field names
% are defined as:
%
% Mandatory Arguments:
%
% 'ImageIn'
%   - the mammogram region
%
% Optional Arguments:
%
% 'Tau'
%   - The expected width of a CLS in pixels. Default 4
%
% 'Len'
%   - The expected ratio of width to length. Default 8
%
% 'NumAngles'
%   - The number of angles at which to apply the filter. Default 12
%
% 'NumScales'
%   - The number of scales at which to apply filters. For each additional
%     scale Tau increases by 2 pixels. Default 1
%
% Return Value:
%
%   MB_CLS_SELECTION returns RESPONSE the magnitude response to the filters
%   and ORIENTATION the orientation that returned maximal response at each
%   pixel
%
% References:
%     'Gabor Filters and Phase Portraits for the Detection of Architectural
%     Distortion in Mammograms' - Rangayyan, R.M.; Ayres, F.J; 2006
%
% Notes: Check the angles of the orientation response closely, it is likely that if
% you are using them with a common co-ordinate system you will want to take
% the negative of the angles. This is because co-ordinate numbering of an
% image in Matlab mean you effectively have an inverted y-axis. e.g. If you
% use quiver to plot the direction of the orientations, the arrows will
% align correctly only if you take the negative first.
%
% Unpack the arguments:
% 
args = u_packargs(varargin, '0', ... % strict mode
		  {...  % The mandatory arguments
          'ImageIn',...
          },... % The optional arguments
          'Tau', 4,...
          'Len', 8,...
          'NumAngles', 12, ...
          'NumScales', 1, ...
          'HighPassFilter', 1, ...
          'Normalise', 1,...
          'Plot', 0,...
          'Debug', false);

clear varargin;

%Get size of image
[rows cols] = size(args.ImageIn);

%Work out angles at which to apply filters
angles = -90:(180/args.NumAngles):90 - (180/args.NumAngles);

%Pre-allocate response/orientations matrices
response = -inf;
response = response(ones(rows, cols));
orientation = zeros(rows, cols);

%Make sure image is double
args.ImageIn = double(args.ImageIn);

for jj = 1:args.NumScales
    
    %make base Gabor filter with zero orientation
    sx = args.Tau / (2*sqrt(2*log(2)));
    sy = args.Len*sx;
    scale = 2*ceil(sy);
    
    %make sure scale doesn't result in a filter bigger than the image
    if 2*scale + 1 >= min([rows cols]);
        warning('Calculated filter size is larger than image. Reducing filter scale'); %#ok
        scale = floor((min([rows cols])-2)/2);
    end
    
    xx = repmat((-scale:scale)', 1 , 2*scale+1);
    yy = repmat(-scale:scale, 2*scale+1, 1);

    gabor_zero = exp(-0.5*((xx.^2 / sx^2) + (yy.^2 / sy^2))).*cos(2*pi*xx / args.Tau)...
        / (2*pi*sx*sy);
    gabor_zero = gabor_zero ./ sum(gabor_zero(:));
      
    if args.HighPassFilter
        %High-pass filter the image
        %make Gaussian filter
        gauss_filt = fspecial('gaussian', scale, sy);
        lowpass_im = imfilter(args.ImageIn, gauss_filt, 'symmetric');
        args.ImageIn = args.ImageIn - lowpass_im;
        if args.Debug
            figure; imagesc(args.ImageIn); axis image; colormap(gray(256));
        end
    end
    
    clear sx sy scale xx yy
    
    %Pre-allocate response/orientations matrices for this scale
    response_scale = -inf;
    response_scale = response_scale(ones(rows, cols));
    orientation_scale = zeros(rows, cols);
    
    

    for ii = 1:args.NumAngles
        
        %make orienatated Gabor filter
        theta = angles(ii);
        gabor_theta = imrotate(gabor_zero, theta, 'bilinear', 'crop');
        %gabor_theta = gabor_theta ./ sum(gabor_theta(:));
        
        %Apply gabor filter to image
        %fr = corrDn(args.ImageIn, gabor_theta, 'reflect1');
        fr = imfilter(args.ImageIn, gabor_theta, 'symmetric');
        
        % Find which pixels have an increased response
        ind = fr > response_scale;
        %display(['changes = ', num2str(sum(ind(:)))]);
        
        % Apply new response and orientation to appropriate pixels
        response_scale(ind) = fr(ind);
        orientation_scale(ind) = angles(ii);
    end
    
    if args.Normalise
    %Normalise the response for the current scale
        response_scale = ((response_scale - min(response_scale(:))) ./ ...
            (max(response_scale(:)) - min(response_scale(:))));
    end
    
    if args.Plot  
        figure; imagesc(response_scale); axis('image'); colormap(gray(256));
    end

    % Apply new response and orientation to appropriate pixels
    inds = response < response_scale;
    response(inds) = response_scale(inds);
    orientation(inds) = orientation_scale(inds);
    
    %Increase Tau for next scale iteration
    args.Tau = args.Tau + 2;
end

%Convert orientation degrees into radians (maybe allow user to choose
%this?)
orientation = orientation * 2*pi / 360;

if args.Plot
    figure; imagesc(response); axis('image'); colormap(gray(256));
    figure; imagesc(orientation); axis('image'); colormap(jet(256));
end

%Main function end