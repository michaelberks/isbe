function [mask, lineprm, nipple] = findNipplePecLine(objLM, im, BS_mask, resize_factor, hough_detsens)
% Based on the image find pectoral line and nipple point
% From Andreas Eilschou, 2010
% im: image
% BS_mask: mask of the whole breast
% resize_factor: if image was resized by a factor, when compared to BS_mask
% hough_detsens - Line detection sensitivity. Between 0 and 1

% Hough transform
[accum, axis_rho, axis_theta, lineprm] = Hough_Grd(im, 8, hough_detsens);

if 0 %houghfigure
    figure; imagesc(axis_theta*(180/pi), axis_rho, accum); axis xy;
    xlabel('Theta (degree)'); ylabel('Rho (pixels)');
    title('Accumulation Array from Hough Transform');
end

mask = +BS_mask;
if numel(lineprm) > 0
    % Upscale the pectoral line to the original size.
    % lineprm contains the line parameters rho,theta.    
    % Adjust rho. Theta is unchanged.
    lineprm(1) = (lineprm(1)-0.5)*resize_factor+0.5;
    
    % Identify the pectoral in the mask
    [~,pectoral_mask,~] = pectoral_area(BS_mask,lineprm(1),lineprm(2));
    mask(pectoral_mask) = objLM.labelPectoral;
       
    % Locate the nipple
    % Move a line parallel to the pectoral line to the right and locate the
    % last breast point. This is faster than computing the point-line distances
    nipple = locateNippleParallelPectoral(mask, lineprm, objLM.labelBreast, objLM.labelBackground);
else
    nipple = [];
end

end

