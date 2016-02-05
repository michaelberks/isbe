function [mask,lineprm,nipple] = fitLeastSquaresLinePectoral(objLM, BS_mask)
% SEGMENTPECTORAL finds the pectoral muscle in the image and returns a mask
% (0 background, 1 breast, 2 pectoral), polar coordinates for the pectoral
% line (rho,theta), and the nipple location.
%
% im            - image to be processed
% BS_mask       - mask specifying air and skin
% hough_detsens - Line detection sensitivity. Between 0 and 1
%                 The lower, the more sensitive.
% Joselene Marques 2014

    % Fit boundary
    pm = BS_mask == objLM.labelPectoral;
    boundaryImg=colfilt(pm,[3 3],'sliding',@(A)any(A==1)&any(A==0)); 
    [y,x]=find(boundaryImg);% boundaryImg(1:50,:) = 0; boundaryImg(:,1:50) = 0;  boundaryImg(end-50:end,:) = 0; boundaryImg(:,end-50:end) = 0;imshow(boundaryImg,[])
    
    % Crop small margin
    % New solution
    unx = unique(x);
    uny = zeros(length(unx),1);
    for i=1:length(unx)
        uny(i) = max(y(x==unx(i)));
    end
    x=unx; y=uny; 
    %----
%     % Old solution
%     idx = ( y > (size(boundaryImg,1)*0.025) & x > (size(boundaryImg,2)*0.025) );
%     x= x(idx);y= y(idx);

    % Fit least square
    coeffsx = polyfit(x, y, 1); 

    % Intersection of the pectoral line perpendicular and the origin
    lp=[coeffsx(1); -1; coeffsx(2)]; % Pectoral line
    xp=cross(cross([lp(1);lp(2);0],[1;1;1]),lp); 
    xp=reshape(xp/xp(3),3,1);
    xprm=xp(1);
    yprm=xp(2);  
    
    % Pectoral line in polar coordinates
    [th,r] = cart2pol(xprm,yprm);
    lineprm = [r th];

    if 0
        figure; imshow(boundaryImg,[]); hold on; 
        x1= min(y); x2= max(y);
        y1 = polyval(coeffsx, x1);
        y2 = polyval(coeffsx, x2);            
        plot(x2, y2,'g*')
        plot(x1, y1,'g*')

        fittedY = polyval(coeffsx, x);
        plot(x, fittedY, 'r*')
        plot(x, y, 'g*')
    end        

    % Identify the pectoral in the mask
    mask = uint8(BS_mask==objLM.labelPectoral | BS_mask==objLM.labelBreast); 
    [~,pectoral_mask,~] = pectoral_area(BS_mask,lineprm(1),lineprm(2));
    mask(pectoral_mask) = objLM.labelPectoral;
    
    nipple = locateNippleParallelPectoral(BS_mask, lineprm, objLM.labelBreast, objLM.labelBackground);

end
