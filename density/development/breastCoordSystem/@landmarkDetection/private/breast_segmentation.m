function [mask,pectoral,nipple] = breast_segmentation(objLM, BS_mask,fig,im_crop,hough_det_sens)
% BREAST_SEGMENTATION finds the pectoral muscle in the image and returns a 
% mask (0 background, 1 breast, 2 pectoral), two points on the pectoral
% line, and the nipple point.
%
% Input:
% BS_mask        - 0 background, 1 breast including pectoral area
% airskin_iter   - Number of iterations in Naga's air-skin segmentation
% hough_det_sens - Sensitivity of line detection in the Hough transform.
%                  The range of the value is from 0 to 1, open ends. The 
%                  smaller the value, the more features in the image will 
%                  be considered as lines.
% fig            - set whether show figure. Number of seconds before next image is 
%                  processed. 
% im_crop        - Number of pixels to crop from borders
%                  [top, bottom, left, right]
%
% Output:
% mask           - 0 background, 1 breast, 2 pectoral
% pectoral       - Two points [I1,J1;I2,J2] on the pectoral line
% nipple         - One point [I1,J1]
%
% Andreas Eilschou, 2010
% Joselene Marques - excluded the automatic determination of direction of
%       the breast (left or right) and the Flip

im = objLM.image;

if nargin < 3
    fig = 0;
end
if nargin < 4
    im_crop = [0 0 0 0];
end
if nargin < 5
    hough_det_sens = 0.08;
end

% Initial smoothing
im = conv2(double(im),fspecial('gaussian'),'same');

% Crop image
im = im(im_crop(1)+1:end-im_crop(2),im_crop(3)+1:end-im_crop(4));
BS_mask = BS_mask(im_crop(1)+1:end-im_crop(2),im_crop(3)+1:end-im_crop(4)); % JM- it was not here... 


% Air-skin segmentation - uses a mean shape
% BS_mask = SaveBSRegion_simplified(im,airskin_iter);

% Pectoral segmentation
%[mask,lineprm,nipple] = segmentPectoral(objLM, im,BS_mask,hough_det_sens);
% Fit least-squares polynom to the pectoral mask
[mask,lineprm,nipple] = fitLeastSquaresLinePectoral(objLM, BS_mask);

if fig > 0
    DrawPolarPectorial(im,  mask, lineprm, nipple)
end

% Convert line parameters (rho,theta) to (x1,y1,x2,y2)
fig = 0;
[mask, pectoral, nipple] = convertPolarNipplePectLine(objLM, mask, lineprm, nipple, im_crop, fig);

% Set background and breast as defined in the Class
mask(mask==0)= -1; % Background;
mask(mask==1) = objLM.labelBreast;
mask(mask==-1) = objLM.labelBackground;

end

