function objCS = setManualLandmarks(objCS, img)
%BRSTPAR pics the breast parameters manually from the image.
%   [XA,YA,XB,YB,XC,YC,XN,YN]=BRSTPAR(IMG) computes the breast parameters 
%   (XA,YA), (XB,YB), (XC,YC), (XN,YN) from the image IMG from for which 
%   reference points are given manually.  
%
%   [...,C1,C2,L]=BRSTPAR(IMG) additionally returns the breast boundary 
%   parabolae C1 and C2 and the pectoral line L. 
    img = processImageMask(img);
    [objCS.xA, objCS.yA ,objCS.xB, objCS.yB, objCS.xC, objCS.yC, ...
        objCS.xM, objCS.yM, objCS.c1n, objCS.c2n, objCS.lp] = brstpar(img);
    objCS.typeCS = 'manual';
end
