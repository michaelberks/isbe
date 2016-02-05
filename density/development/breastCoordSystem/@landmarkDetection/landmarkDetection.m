classdef landmarkDetection 

    properties (GetAccess = public, SetAccess = public)
        subjID = '';
        mask = [];
        pectoral = [];
        nipple = [];
        pctLine = [];
        image = [];
        
        labelPectoral = 2;
        labelBreast = 1;    
        labelBackground = 3;        
        flipped = 0;
    end

    
  % ----------------------------  Public  ---------------------------------  
  methods (Access = public) 
    
      function objLM = setLabels(objLM,labelBackground, labelBreast, labelPectoral)       
        objLM.labelBackground = labelBackground; 
        objLM.labelBreast = labelBreast; 
        objLM.labelPectoral = labelPectoral; 
      end
    
    % Class constructor
    function objLM = landmarkDetection(I, M,labelBackground, labelBreast, labelPectoral)
%         filepath = which('BreastCoodSystem');
%         idx = strfind(filepath, '@');
%         addpath([filepath(1:idx-1) 'matfiles'])
           
        objLM.image = processImageMask(I);
        objLM.mask = processImageMask(M);
        
        objLM.labelBackground = labelBackground; 
        objLM.labelBreast = labelBreast; 
        objLM.labelPectoral = labelPectoral; 
        
            
        % Flip to Left if it is necessary    
        [objLM.image, objLM.mask, objLM.flipped ] = detectMaskDirecFlipImage2Left(objLM.image, objLM.mask, objLM.labelBackground);
    end
      
    
    function objLM = detectCSLandmarks(objLM, bFig)

        % Resize and ensure image and mask are the same size
        tmpimage = imresize(objLM.image,size(objLM.mask,1)/size(objLM.image,1)); 
        objLM.image = zeros(size(objLM.mask)); 
        objLM.image(:) = tmpimage(1:length(objLM.image(:)));
        
        % Segmentation
        if ~exist('bFig','var')
            bFig = 0;
        end
        
        % Crops mask to ensure parabolae will not be misguided
        % Target code: matrix(im_crop(1)+1:end-im_crop(2),im_crop(3)+1:end-im_crop(4))
        if ~exist('im_crop','var')
            im_crop = [0 0 0 0]; 
        end
        [objLM.mask,objLM.pectoral,objLM.nipple] = breast_segmentation(objLM, objLM.mask,...
            bFig, im_crop);

        objLM.pctLine = LineConv(objLM.pectoral); 
    end    
   
    % bTypeFigure: 1=image, 2=mask, 3=both
    function showSegmentation(objLM, bTypeFigure)
        if bTypeFigure==3, subplot(1,2,1); end
        if (bTypeFigure==3) || (bTypeFigure==2)
            imshow(objLM.mask,[]); 
        end
        if bTypeFigure==3, subplot(1,2,2); end
        if (bTypeFigure==3) || (bTypeFigure==1)
            imshow(objLM.image,[]);                                    
        end
        hold on; 
        showPoints(objLM)
    end
        
    function showPoints(objLM)
        plot(objLM.pectoral(1,2),objLM.pectoral(1,1),'r*');
        plot(objLM.pectoral(2,2),objLM.pectoral(2,1),'r*');
        plot(objLM.nipple(2),objLM.nipple(1),'g*');                          
    end

  
  end

  % --------------------------------------------------------------------------------------------
  methods (Access = public, Static)
      
    % Based on the image find pectoral line and nipple point
    % im: image
    % BS_mask: mask with breast and pectoral area 
    % hough_detsens - Line detection sensitivity. Between 0 and 1  
%     function [pctLine, nipple] = findNipplePectoralLine(objLM, im, hough_detsens, fig)
%         tmpMask = (im==1 | im==2);
%         im(im==2) = 255;
%         im(im==1) = 2;
%         
%         [tmpMask, lineprm, nipple] = findNipplePecLine(objLM, im, tmpMask, 1, hough_detsens);
%         if fig > 0
%             DrawPolarPectorial(im,  tmpMask, lineprm, nipple)
%         end
%         
%         % Convert line parameters (rho,theta) to (x1,y1,x2,y2)
%         [~, tmpPectoral, nipple] = convertPolarNipplePectLine(tmpMask, lineprm, nipple, [0 0 0 0 ], fig);
% 
%         pctLine = LineConv(tmpPectoral); 
%     end
  end
  
end
