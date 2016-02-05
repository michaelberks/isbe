function [mask, bFlipped] = detectFlipMask2Left(mask, labelBackground)
% Automatically determine the direction of the breast (left or right)
% and flip the right images to left
    
    im_direction =  maskDirection(mask,labelBackground);

    % Flip image
    bFlipped = strcmp('right',im_direction);
    if bFlipped
        mask = fliplr(mask);
    end
end