function [corrC, pCorr, X1, Y1, X2, Y2] = evaluateCoodSystem(objCS1, objCS2, im1, im2, mask1, corrMethod, bShowFigure)
% returns the correlation coefficient and the pValue

    % Load dicom 
    im1 = processImageMask(im1);
    im2 = processImageMask(im2);
    
    im1Size = size(im1);
    im2Size = size(im2);
    

    % Random points from image 1
    if length(unique(mask1)) > 2
        points = find(mask1 == objCS1.labelBreast); 
    else
        points = find(mask1); 
    end
    idx = randi(length(points), [1, objCS1.nrPtsEvaluate]);
    points = points(idx);    
    [X1, Y1] = ind2sub(im1Size,points);

    % Coorrespondent points in Breast Coordinate System
    [S,PHI] = objCS1.cart2CoordSystem(X1, Y1);
    
    % Cartesian correspondent points in Image2
    [X2, Y2] = objCS2.coordSystem2Cart(S, PHI);

    % Exclude points without correspondence, for example, correspondent is 
    % outside the image
    [X1, Y1, X2, Y2] = exclCoodSystemOutImage(im2Size, X1, Y1, round(X2), round(Y2));
    
    intensities1 = single(im1(sub2ind(im1Size,X1,Y1)));
    intensities2 = single(im2(sub2ind(im2Size,X2,Y2)));
    
    if ~isempty(intensities1)
        [corrC,pCorr] = corr(intensities1,intensities2,'Type',corrMethod);
    else
        corrC = NaN;
        pCorr = NaN;
    end
    
    % Plot image and points
    if bShowFigure
        figure; 
        subplot(1,2,1); imshow(im1,[]); hold on; plot(Y1, X1, 'g*'); 
        subplot(1,2,2); imshow(im2,[]); hold on; plot(Y2, X2, 'g*')
    end
end


function [X1, Y1, X2, Y2] = exclCoodSystemOutImage(imSize, X1, Y1, X2, Y2)
    idx =  isnan(X2) | isnan(Y2) | Y2<1 | X2<1 | X2>imSize(1) | Y2>imSize(2) | imag(X2(:))~=0 | imag(Y2(:))~=0;    
    X1 = X1(~idx);
    X2 = X2(~idx);
    Y1 = Y1(~idx);
    Y2 = Y2(~idx);
end



