function [I1, M1, objCS1] = getCoordSystem(objData,  I1, M1, bFig)      
    if ~exist('bFig','var')
        bFig = 0;
    end
    
    % right images have to be flipped
    [M1, flipped] = detectFlipMask2Left(M1, objData.labelBackground);
    if flipped
        I1 = fliplr(I1);
    end
    objLM = landmarkDetection(I1, M1,objData.labelBackground, objData.labelBreast, objData.labelPectoral);
    objLM = objLM.detectCSLandmarks(bFig); 
    objCS1 = BreastCoodSystem();   
    objCS1 = objCS1.setLabels(objData.labelBackground, objData.labelBreast, objData.labelPectoral); 
    objCS1 = objCS1.mask2BreastCoodSystem(objLM.mask, objLM.nipple, objLM.pctLine);
    objCS1.flipped = flipped;
end