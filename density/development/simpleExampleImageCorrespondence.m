function simpleExampleImageCorrespondence
    % Example of how to use the imageCorrespondentPoints funtion. This function gets the masks and
    % images of both time points and  the sample points from one of the images. The output
    % is the correspondent points in the other image. It is possible to give as input the constant
    % value used as breast, background and pectoral areas, as well.
    % BCS
    % Input data: the masks for the images containing the breast and pectoral area segmentation. 
    %   Additionally, the constant value used as breast, background and pectoral areas in the masks 
    %   should also be informed.
    % Flipping: The BCS only computes left-sided images, so if the breast is in the right side
    %   (normally RMLO), it is flipped automatically. The correct image and mask (flipped or not) is
    %   available in the class Landmarks. 
    % This version is only applicable for MLO views.
    % Author : Joselene Marques

    % Add path to Coord System code
    addpath('breastCoordSystem');     
    
    % read data
    dcInfo{1} = dicominfo('image1'); % ASSURE_MANC_50-150TrainingSet_Cancer-00008-RMLO150347-20121204
    dcInfo{2} = dicominfo('image2'); % ASSURE_MANC_50-150TrainingSet_Cancer-00009-RMLO100419-20120912        
    resLoad = load('masksManc','masks');  
    masks = resLoad.masks;
    
    
    % Random points from image 1 
    masktmp = fliplr(masks{1}); % Masks and images will be flipped to Left in the BCS code
    labelBreast = 3;
    nrPtsEvaluate = 150;
    if length(unique(masks{1})) > 2
        Points_mask1 = find(masktmp == labelBreast); 
    else
        Points_mask1 = find(masktmp); 
    end
    idx = randi(length(Points_mask1), [1, nrPtsEvaluate]);
    Points_mask1 = Points_mask1(idx);    
    
    %--------------------------------------------------------------------------------
    %        The important part - This line you add in your code ! 
    % -----------------------------------------------------------------------------    
    % In this example, masks and images will be flipped to Left, since they are Rmlo!
    bShowCoordSystemFigures = 1;
    [bcs1, bcs2] = BreastCoodSystem.imageCorrespondentPoints(dcInfo{1}, dcInfo{2}, masks{1}, masks{2}, Points_mask1, bShowCoordSystemFigures);

    
    %--------------------------  Just an example of application 
    im1 = fliplr(processImageMask(dcInfo{1}));
    im2 = fliplr(processImageMask(dcInfo{2}));
    intensities1 = single(im1(sub2ind(size(masks{1}),bcs1.X,bcs1.Y)));
    intensities2 = single(im2(sub2ind(size(masks{2}),bcs2.X,bcs2.Y)));    
    if ~isempty(intensities1)
        [corrC,pCorr] = corr(intensities1,intensities2);
    else
        corrC = NaN;
        pCorr = NaN;
    end
    fprintf('Correlation and p-Value: %1.2f  (%1.4f) \n', corrC, pCorr);
    % Plot image and points
    figure; 
    subplot(1,2,1); imshow(im1,[]); hold on; plot(bcs1.Y, bcs1.X, 'g*'); 
    subplot(1,2,2); imshow(im2,[]); hold on; plot(bcs2.Y, bcs2.X, 'g*')

end
