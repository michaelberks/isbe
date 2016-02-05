function simpleExampleAutoCSdicom
    % Example of how to use the automatic version of the Breast Coordinate System.
    % BCS
    % Input data: For computation, one needs the masks for the images containing the breast and
    %   pectoral area segmentation. Additionally, the constant value used as breast, background 
    %   and pectoral aresa in the masks should be informed.
    % Flipping: The BCS only computes left-sided images, so if the breast is in the right side
    %   (normally RMLO), it is flipped automatically. The correct image and mask (flipped or not) is
    %   available in the class Landmarks. 
    % This version is still only applicable for MLO views.
    % Author : Joselene Marques
    
    set(0,'DefaultFigureWindowStyle','docked')
    warning('off', 'images:initSize:adjustingMag');
    warning('off', 'images:imshow:magnificationMustBeFitForDockedFigure');

    % Mask's labels 
    labelPectoral = 2;
    labelBreast = 3;    
    labelBackground = 1;    
    bFig = 1; %%%new - Exibe figuras durante Coord. System detection 
    
    % read data
    dcInfo{1} = dicominfo('image1'); % ASSURE_MANC_50-150TrainingSet_Cancer-00008-RMLO150347-20121204
    dcInfo{2} = dicominfo('image2'); % ASSURE_MANC_50-150TrainingSet_Cancer-00009-RMLO100419-20120912        
    resLoad = load('masksManc','masks'); %%%new  
    masks = resLoad.masks; %%%new
    
    % Add path to Coord System code
    addpath('breastCoordSystem');     
    
    % Example with manual annotation
    fprintf('\n\nExample with automatic segmentation \n')
    objCS{length(dcInfo)} = [];
    objLM{length(dcInfo)} = []; %%%new
    for s=1:length(dcInfo)        
        objLM{s} = landmarkDetection(dcInfo{s}, masks{s},labelBackground, labelBreast, labelPectoral); %%%new
        objLM{s} = objLM{s}.detectCSLandmarks(bFig); %%%new
        objCS{s} = BreastCoodSystem();  
        objCS{s} = objCS{s}.setLabels(labelBackground, labelBreast, labelPectoral); 
        objCS{s} = objCS{s}.mask2BreastCoodSystem(objLM{s}.mask, objLM{s}.nipple, objLM{s}.pctLine);    
    
        figure;
        % The automatic Coord. System requires image on left side, so the original image was flipped
        objCS{s}.showCoodSystem(objLM{s}.image); %%%new - Image from ladmark object, since it might be flipped
    end    
    
    % Example of mapping: evaluate the correlation coefficient of intensities
    % between two breast images        
    objCS{1}.nrPtsEvaluate = 1000; % Number of evaluation's points 
    [corrC, pCorr] = evaluateCoodSystem(objCS{1}, objCS{2}, objLM{1}.image, objLM{2}.image, objLM{1}.mask, 'Pearson', 1);
    fprintf('Correlation and p-Value: %1.2f  (%1.4f) \n', corrC, pCorr);
    
end
