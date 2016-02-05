function [bcs1, bcs2] = imageCorrespondentPoints(dicomIm1, dicomIm2, mask1, mask2, Points_mask1, ...
    bShowFigure, labelBackground, labelBreast, labelPectoral )
    % Points_mask1: logical matrix with same size as mask1 and dicomIm1

    if ~exist('bShowFigure','var')
        bShowFigure = 0;
    end
    if ~exist('labelBackground','var')
        labelBackground = 1;
    end
    if ~exist('labelBreast','var')
        labelBreast = 3;
    end
    if ~exist('labelPectoral','var')
        labelPectoral = 2;
    end
    
    
    % Image 1
    objLM = landmarkDetection(dicomIm1, mask1,labelBackground, labelBreast, labelPectoral); 
    objLM = objLM.detectCSLandmarks(bShowFigure); 
    objCS1 = BreastCoodSystem();  
    objCS1 = objCS1.setLabels(labelBackground, labelBreast, labelPectoral); 
    objCS1 = objCS1.mask2BreastCoodSystem(objLM.mask, objLM.nipple, objLM.pctLine);   
    bcs1.flipped = objLM.flipped | objCS1.flipped;
    %figure;
    % The automatic Coord. System requires image on left side, so the original image was flipped
    if bShowFigure
        objCS1.showCoodSystem(objLM.image);
    end
    
    % Image 2
    objLM = landmarkDetection(dicomIm2, mask2,labelBackground, labelBreast, labelPectoral); 
    objLM = objLM.detectCSLandmarks(bShowFigure); 
    objCS2 = BreastCoodSystem();  
    objCS2 = objCS2.setLabels(labelBackground, labelBreast, labelPectoral); 
    objCS2 = objCS2.mask2BreastCoodSystem(objLM.mask, objLM.nipple, objLM.pctLine); 
    bcs2.flipped = objLM.flipped | objCS2.flipped;
    if bShowFigure
        objCS2.showCoodSystem(objLM.image);
    end
    
    % Load dicom 
    dicomIm1 = processImageMask(dicomIm1);
    dicomIm2 = processImageMask(dicomIm2);
    
    %im1Size = size(dicomIm1);
    im2Size = size(dicomIm2);
      
    [X1, Y1] = ind2sub(size(mask1),Points_mask1);

    % Coorrespondent points in Breast Coordinate System
    [S,PHI] = objCS1.cart2CoordSystem(X1, Y1);
    
    % Cartesian correspondent points in Image2
    [X2, Y2] = objCS2.coordSystem2Cart(S, PHI);

    % Exclude points without correspondence, for example, correspondent is 
    % outside the image
    [X1, Y1, X2, Y2] = exclCoodSystemOutImage(im2Size, X1, Y1, round(X2), round(Y2));
    
    bcs1.X = X1;
    bcs1.Y = Y1;
    bcs2.X = X2;
    bcs2.Y = Y2;    
    
    % Plot image and points
    if bShowFigure
        figure; 
        subplot(1,2,1); imshow(dicomIm1,[]); hold on; plot(Y1, X1, 'g*'); 
        subplot(1,2,2); imshow(dicomIm2,[]); hold on; plot(Y2, X2, 'g*')
    end
end


function [X1, Y1, X2, Y2, idx] = exclCoodSystemOutImage(imSize, X1, Y1, X2, Y2)
    idx =  isnan(X2) | isnan(Y2) | Y2<1 | X2<1 | X2>imSize(1) | Y2>imSize(2) | imag(X2(:))~=0 | imag(Y2(:))~=0;    
    X1 = X1(~idx);
    X2 = X2(~idx);
    Y1 = Y1(~idx);
    Y2 = Y2(~idx);
end



