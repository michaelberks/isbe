function [mask,lineprm,nipple] = segmentPectoral(objLM, im,BS_mask,hough_detsens)
% SEGMENTPECTORAL finds the pectoral muscle in the image and returns a mask
% (0 background, 1 breast, 2 pectoral), polar coordinates for the pectoral
% line (rho,theta), and the nipple location.
%
% im            - image to be processed
% BS_mask       - mask specifying air and skin
% hough_detsens - Line detection sensitivity. Between 0 and 1
%                 The lower, the more sensitive.
%
% Andreas Eilschou, 2010


%%%% PARAMETERS %%%%
resize_height = 200;

% resize images for faster computation
resize_im = imresize(im,[resize_height NaN]);
resize_BS_mask = double(imresize(BS_mask,[resize_height NaN]));
resize_factor = size(im,1)/resize_height;

% % the center of gravity (CG)
% STATS = regionprops(resize_BS_mask, 'Centroid');
% CG = round(STATS(1).Centroid);

BS_mask = (BS_mask==objLM.labelPectoral | BS_mask==objLM.labelBreast);
% Find pectoral line inside the pectoral area in the mask - JM
jmAdjustingPectoral = 0;
if jmAdjustingPectoral
    resize_BS_mask = resize_BS_mask==objLM.labelPectoral;

    % the center of gravity (CG)
    STATS = regionprops(resize_BS_mask, 'Centroid');
    CG = round(STATS(1).Centroid);

    se = strel('square',30);
    resize_BS_mask = imdilate(resize_BS_mask,se);
end


% reduce the ROI to above a line through CG with slope 1.5.
% CG(1) is horizontal, CG(2) is vertical
if ~jmAdjustingPectoral
    slope = 1.5; 
    [N,M] = size(resize_BS_mask);
    for j = 1:M
        line = CG(2)+floor((CG(1)-j)*slope);
        if (line < 1)
            line = 1;
        elseif (line > N)
            line = N;
        end
        for i = 1:line % from  line:N 
            resize_BS_mask(i,j) = 0;
        end
    end
end

% Find edges using the 3x3 Sobel operator
SOBEL = [1 2 1; 0 0 0; -1 -2 -1];
Gy = conv2(double(resize_im),SOBEL,'same');
Gx = conv2(double(resize_im),SOBEL','same');

% Gradient magnitude
G = sqrt(Gx.^2+Gy.^2) .* resize_BS_mask;

% Based on the image find pectoral line and nipple point
[mask, lineprm, nipple] = findNipplePecLine(objLM, G, BS_mask, resize_factor, hough_detsens);


end
