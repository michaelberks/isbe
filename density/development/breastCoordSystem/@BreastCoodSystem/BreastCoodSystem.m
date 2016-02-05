classdef BreastCoodSystem 

    properties (GetAccess = public, SetAccess = private)
        typeCS = '';
    end
    
    properties (GetAccess = public, SetAccess = public)
        subjID = '';
        xA = [];
        yA = [];
        xB = [];
        yB = [];
        xC = [];
        yC = [];
        xM = [];
        yM = [];
        c1n = [];
        c2n = [];
        Cpoints = [];
        lp = [];
        
        % How many lines to show in the coord system figure
        resolution = 12; %36
        
        labelPectoral = 2;
        labelBreast = 1;    
        labelBackground = 3; 
        nrPtsEvaluate = 1000;
        flipped = 0; % Image and mask were flipped.

    end
      
    
  % ----------------------------  Public  ---------------------------------  
  methods (Access = public) 
    

    function objCS = setLabels(objCS,labelBackground, labelBreast, labelPectoral)       
        objCS.labelBackground = labelBackground; 
        objCS.labelBreast = labelBreast; 
        objCS.labelPectoral = labelPectoral; 
    end
      
    % Class constructor:
    function objCS = BreastCoodSystem()

    end
    
    % Generate the breast coordinate system.
    % INPUT:
    % mask: breast mask where 1=breast.
    % nipple: coordenates x and y of the nipple point.
    % pctLine: pectorial line.   
    % OUTPUT:
    % objCS: BreastCoodSystem object.
    function objCS = mask2BreastCoodSystem(objCS, mask, nipple, pctLine)
        
        %objCS.mask = mask;
        
        % Flip to Left if necessary
        [mask, objCS.flipped] = detectFlipMask2Left(mask, objCS.labelBackground);  
        
        % Fuction to convert Nipple point, pectoral line into
        % BC cordinate landmarks         
        breastCS = mask2par(objCS, mask, nipple(2), nipple(1), pctLine'); 
    
        objCS.xA = breastCS.xA;
        objCS.yA = breastCS.yA;
        objCS.xB = breastCS.xB;
        objCS.yB = breastCS.yB;
        objCS.xC = breastCS.xC; 
        objCS.yC = breastCS.yC;        
        objCS.c1n = breastCS.c1n;
        objCS.c2n = breastCS.c2n;
        objCS.lp = pctLine';
        objCS.Cpoints = breastCS.Cpoints;
        objCS.typeCS = 'auto';
    end
        
    showCoodSystem(objCS, image)
    
    showPoints(objCS)
    
    % function: coordSystem2Cart
    % convert points from breast coordinate system to cartesian.
    % INPUT:
    % objCS: BreastCoodSystem object.
    % S: relative distance(s) from the nipple along the parabolic line
    % (array or scalar).
    % PHI: the direction(s) of the parabola at the nipple (array or
    % scalar).
    % OUTPUT:
    % X, Y: coordenates in cartesian system.
    % J: the Jacobian of the coordinate transform.
    function [X, Y, J] = coordSystem2Cart(objCS, S, PHI)       
         [Y,X,J] = brst2cart(S,PHI,objCS.xA,objCS.yA,objCS.xB, ...
             objCS.yB,objCS.xC,objCS.yC,objCS.c1n,objCS.c2n); 
    end
    
    % function: coordSystem2Cart
    % convert points from breast coordinate system to cartesian.
    % INPUT:
    % objCS: BreastCoodSystem object.
    % X, Y: coordenates in cartesian system.    
    % OUTPUT:    
    % S: relative distance(s) from the nipple along the parabolic line
    % (array or scalar).
    % PHI: the direction(s) of the parabola at the nipple (array or
    % scalar).
    function [S,PHI] = cart2CoordSystem(objCS, X, Y)       
        [S,PHI] = cart2brst( cast(Y,'double') , cast(X,'double'), objCS.xA,objCS.yA,objCS.xB,objCS.yB,...
             objCS.xC,objCS.yC,objCS.c1n,objCS.c2n); 
    end
    
    % Correlation of intensities
    %[corrC,pCorr] = evaluateCoodSystem(objCS1, objCS2, im1, im2, gridResolution, corrMethod)
    
    function plotConics(objCS)
        hold on;
        plotconic(objCS.c1n);
        plotconic(objCS.c2n);
    end
        
  end
  
  
  % ----------------------------------------------------------------------------
  methods (Access = public, Static)
      [X2, Y2, X1, Y1] = imageCorrespondentPoints(dicomIm1, dicomIm2, mask1, mask2, Points_mask1, ...
           bShowFigure, labelBackground, labelBreast, labelPectoral )

  end
  
  
end
