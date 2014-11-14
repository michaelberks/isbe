function pixel_order = mb_compute_pixel_order(filled_image, synthesis_mode, patch_size)
%
% MB_COMPUTE_PIXEL_ORDER Synthesise an image using an ergodic GMM.
%
% Input:
%  'filled_image'
%  - 
%
% 'synthesis_mode'
% -
%
% Return Value:
%
%   pixel_order
%
% See also:
%
%pre-allocate pxiel_order

if nargin < 2
    synthesis_mode = 'pixel';
end

pixels_to_fill = sum(~filled_image(:));

switch synthesis_mode
    case 'pixel'
        
        start_pix = 1;
        pixel_order = zeros(pixels_to_fill, 1);
        while start_pix <= pixels_to_fill;
            
            %Get the perimeter pixels for this onion-skin
            unfilled_list = get_unfilled_perimeter_pixels(filled_image);
            
            %Randomly permute the pixel order and add to list
            end_pix = start_pix + length(unfilled_list) - 1;
            pixel_order(start_pix:end_pix) = unfilled_list(randperm(length(unfilled_list)));
            
            %Update filled_image and pixel counter
            filled_image(unfilled_list) = 1;
            start_pix = end_pix + 1;
        end
        
    case 'patch'
        
        overlap = (patch_size - 1) / 2;
        pixel_order = [];
        while pixels_to_fill;
            
            %Get unfilled perimeter pixels
            unfilled_list = get_unfilled_perimeter_pixels(filled_image);
            
            %Choose one at random and add to pixel order
            if length(unfilled_list) > 1
                idx = randsample(unfilled_list, 1);
            else
                idx = unfilled_list;
            end
            pixel_order(end+1, 1) = idx; %#ok
            
            %Update filled_image and pixel counters
            [row col] = ind2sub(size(filled_image), idx);
            filled_image(row-overlap:row+overlap,col-overlap:col+overlap) = 1;
            pixels_to_fill = sum(~filled_image(:));
        end
        
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [unfilled_list] = get_unfilled_perimeter_pixels(filled_image)
%
% Get a list of unfilled pixels neighbouring the pixels already filled.
% This is returned as a list of indices -- use ind2sub to get the image
% locations back.

%unfilled_image = bwmorph(filled_image, 'dilate'); % dilation with SE=ones(3)
if any(filled_image(:))
    unfilled_image = imdilate(filled_image, strel('disk', 1));
    unfilled_image = double(unfilled_image) - double(filled_image);
    unfilled_list = find(unfilled_image==1);
else
    unfilled_image = filled_image;
    unfilled_image([1 end], [1 end]) = 1;
    unfilled_list = find(unfilled_image);
end

