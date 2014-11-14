function [CLS_result] = mb_cls_selection(varargin)
%
% MB_CLS_SELECTION Produce binary map of CLS structures from a mammographic
% image
%
% This function selects the CLS structures present in a mammogram region
% using the Gabor filtering method described by Rangayyan and Ayres 2006.
%
% MB_CLS_SELECTION uses the U_PACKARGS interface function
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
% 'GaborFilterArgs'
%   - structure of arguments to pass to the Gabor filter. If left empty
%   Gabor filter function will use its default arguments
%
% 'IgnoreMap'
%   - Binary map, specifying regions to ignore in the image (e.g. off the
%   breast edge or the pectoral muscle
%
% 'GaussDeriv'
%   - 2d array specifying a first derivative of Gaussian filter to find the
%   image gradient. Assumed to be triangular to be applied in x-direction
%   then transpose applied in y-direction
%
% 'SaveFile' - If non-empty, the pathname to save the returned structured
% to
%
% Return Value:
%
%   MB_CLS_SELECTION returns CLS_RESULT a structure containing the
%   following fields:
%    'CLS' 
%      - Binary map of same size as ImageIn, containing 1s at
%        CLS pixels and 0s otherwise
%
%    'GaborResponse' 
%      - magnitude response returned Gabor filtering
%
%    'GaborOrientation' 
%      - orientations returned from Gabor filtering
%
%    'NMS' 
%      - the non-maximally suppressed response to the gabor filter
%
%    'DiscardMap'
%      - the pixels ignored after comparison to the image gradient (NB.
%      Also contains any region specified in IgnoreMap
%
% References:
%     'Gabor Filters and Phase Portraits for the Detection of Architectural
%     Distortion in Mammograms' - Rangayyan, R.M.; Ayres, F.J; 2006

% Unpack the arguments:
args = u_packargs(varargin, '0', ... % strict mode
		  {...  % The mandatory arguments
          'ImageIn',...
          },... % The optional arguments
          'GaborFilterArgs', [],...
          'IgnoreMap', [],...
          'GaussDeriv', [], ...
          'GradientThreshold', 0,...
          'GradientAlignment', 0,...
          'GaborThreshold', 0,...
          'NMS', true,...
          'Thin', 0,...
          'Connectivity', 0,...
          'MinLength', 12, ...
          'FullResult', true,...
          'Plot', 0,...
          'SaveFile', [],...
          'Debug', false);

clear varargin;

if isempty(args.IgnoreMap)
    args.IgnoreMap = false(size(args.ImageIn));
end

%Gabor filter the input image
args.GaborFilterArgs.ImageIn = args.ImageIn;
[CLS_result.GaborResponse CLS_result.GaborOrientation] = ...
    gabor_filter(args.GaborFilterArgs);

CLS_result.GaborOrientation = -CLS_result.GaborOrientation;

%Apply gradient threshold if selected
if args.GradientThreshold || args.GradientAlignment
    
    %Build Gaussian Derivative filter if not supplied
    if isempty(args.GaussDeriv)
        %Make Gaussian filter
        sigma = 5;
        width = 21; %Should probably give the user the option of choosing this

        ssq = sigma^2;
        [x,y] = meshgrid(-width:width,-width:width);
        args.GaussDeriv = -x.*exp(-(x.*x+y.*y)/(2*ssq))/(pi*ssq);
    end
    %Calulate image gradient and normal direction to gradient
    grad_x = imfilter(args.ImageIn, args.GaussDeriv, 'conv','replicate');
    grad_y = imfilter(args.ImageIn, args.GaussDeriv', 'conv','replicate');
end

%%%
if args.Debug
    temp = args.ImageIn;
    temp(args.IgnoreMap) = 0;
    figure; imagesc(temp); axis image; clear temp;
    db_axes(1) = gca;
    c1 = caxis;
end
%%%

if args.NMS
    gabor_nms = CLS_result.GaborResponse;
    gabor_nms(args.IgnoreMap) = 0;
    
    % Apply non-maximal suppression
    CLS_result.NMS = ...
        mb_non_maximal_supp(gabor_nms, CLS_result.GaborOrientation);
    
    if args.FullResult
        CLS_result.NMSMap = ...
            ~args.IgnoreMap & ~CLS_result.NMS;
    end
    
    args.IgnoreMap = args.IgnoreMap | ~CLS_result.NMS;
    
    %%%
    if args.Debug
        temp = args.ImageIn;
        temp(args.IgnoreMap) = 0;
        figure; imagesc(temp); axis image; clear temp; caxis(c1);
        db_axes(end+1) = gca; %#ok
    end
    %%%
end

if args.GradientAlignment
    % Rangayyan: compute gradient normals and discard based on orientation
    normal_to_gradient = atan(-grad_x ./ grad_y);
    grad_diff = mod(normal_to_gradient - CLS_result.GaborOrientation, pi);
    
    args.IgnoreMap = args.IgnoreMap | ...
         (grad_diff < args.GradientAlignment) | (grad_diff > (pi-args.GradientAlignment));
    
    if args.FullResult
        CLS_result.GradientAlignmentMap = ...
            (grad_diff < args.GradientAlignment) | (grad_diff > (pi-args.GradientAlignment));
    end
    
    %%%
    if args.Debug
        temp = args.ImageIn;
        temp(args.IgnoreMap) = 0;
        figure; imagesc(temp); axis image; clear temp; caxis(c1);
        db_axes(end+1) = gca; %#ok
    end
    %%%
end

if args.GradientThreshold
    %discard based on percentage threshold of magnitude
    grad_mag = grad_x.^2 + grad_y.^2;
    args.IgnoreMap = args.IgnoreMap | ...
        (grad_mag > (args.GradientThreshold * max(grad_mag(:))));
    if args.FullResult
        CLS_result.GradientThresholdMap = ...
            (grad_mag > (args.GradientThreshold * max(grad_mag(:))));
    end
    
    %%%
    if args.Debug
        temp = args.ImageIn;
        temp(args.IgnoreMap) = 0;
        figure; imagesc(temp); axis image; clear temp; caxis(c1);
        db_axes(end+1) = gca; %#ok
    end
    %%%
end

if args.GaborThreshold
    % Discard pixels from the Gabor filter response and scale the resulting
    % output between 0 and 1
    gabor_thresh = CLS_result.GaborResponse;
    gabor_thresh(args.IgnoreMap) = 0;
    gabor_thresh = (gabor_thresh - min(gabor_thresh(:))) / max(gabor_thresh(:));    

    % Threshold the output
    args.IgnoreMap = args.IgnoreMap | ...
    	(gabor_thresh < args.GaborThreshold);
    if args.FullResult
        CLS_result.GaborThresholdMap = ...
            (gabor_thresh < args.GaborThreshold);
    end
    clear gabor_thresh;
    
    %%%
    if args.Debug
        temp = args.ImageIn;
        temp(args.IgnoreMap) = 0;
        figure; imagesc(temp); axis image; clear temp; caxis(c1);
        db_axes(end+1) = gca; %#ok
    end
    %%%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Previous steps operated on filter response map, now perform binary
% operations

%Current binary map of CLS pixels is inverse of those we've discarded
cls_map = ~args.IgnoreMap;

if args.Thin
    %Thin the image and remove any 3+ connected pixels
    new_map = bwmorph(cls_map, 'thin', Inf);
    if args.FullResult
        CLS_result.ThinMap = ...
            cls_map & ~new_map;
    end
    cls_map = new_map; clear new_map;
    
    %%%
    if args.Debug
        temp = args.ImageIn;
        temp(~cls_map) = 0;
        figure; imagesc(temp); axis image; clear temp; caxis(c1);
        db_axes(end+1) = gca; %#ok
    end
    %%%
end

% Remove pixels that are 8-connected to more than 2 pixels
if args.Connectivity
    
    temp = cls_map;
    temp(:,[1 end]) = 0; temp([1 end], :) = 0;
    [r_pts c_pts] = find(temp); %find all non-zero and non-boundary pixels in nms
    clear temp;
    
    if args.FullResult
        CLS_result.ConnectivityMap = false(size(cls_map));
    end
    
    for idx = 1:length(r_pts);
        r = r_pts(idx);
        c = c_pts(idx);
        local_win = cls_map(r-1:r+1,c-1:c+1);
        if sum(local_win(:)) > 3
            cls_map(r,c) = 0;
            if args.FullResult
                CLS_result.ConnectivityMap(r,c) = true;
            end
        end
    end
    
    %%%
    if args.Debug
        temp = args.ImageIn;
        temp(~cls_map) = 0;
        figure; imagesc(temp); axis image; clear temp; caxis(c1);
        db_axes(end+1) = gca; %#ok
    end
    %%%
end

if args.MinLength
    %Label image and remove objects of fewer than 10 connected pixels
    cls_label = bwlabel(cls_map, 8);
    
    if args.FullResult
        CLS_result.MinLengthMap = false(size(cls_map));
    end
    
    for i = 1:max(cls_label(:));
        ind = find(cls_label == i);
        if length(ind) < args.MinLength
            cls_map(ind) = 0;
            if args.FullResult
                CLS_result.MinLengthMap(ind) = true;
            end
        end
    end
    
    %%%
    if args.Debug
        temp = args.ImageIn;
        temp(~cls_map) = 0;
        figure; imagesc(temp); axis image; clear temp; caxis(c1);
        db_axes(end+1) = gca; %#ok
    end
    %%%
end

%%%
if args.Debug
    linkaxes(db_axes)
end
%%%

CLS_result.CLS = cls_map;
CLS_result.DiscardMap = args.IgnoreMap;

if ~isempty(args.SaveFile)
    save(args.SaveFile, 'CLS_result');
end

if args.Plot
    figure; imagesc(nms); axis image; colormap(gray(256));
    figure; imagesc(args.ImageIn);
    colormap(gray(256)); axis image; hold on;
    [r c] = find(CLS_result.CLS);
    plot(c,r,'r.', 'MarkerSize', 4);
end