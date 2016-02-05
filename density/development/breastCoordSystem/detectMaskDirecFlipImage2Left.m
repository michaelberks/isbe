function [image, mask, bFlipped] = detectMaskDirecFlipImage2Left(image, mask, labelBackground )
% Automatically determine the direction of the breast (left or right)
% and flip the right images to left

    [mask, bFlipped] = detectFlipMask2Left(mask, labelBackground);
    if bFlipped
        image = fliplr(image);
    end

end