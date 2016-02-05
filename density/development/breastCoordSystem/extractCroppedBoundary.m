function [ y, x, mask] = extractCroppedBoundary(mask,labelBreast,labelBackground)
% Crop the mask to inside a circle that circumscribes the pectoral line and returns the circumscribed boundary
% extractCroppedBoundary(mask,objCS.labelBreast, objCS.labelBackground,xp(2),xp(1))

    if 0 % Old version
        mask = uint8(mask);
        boundaryImg=colfilt(mask,[3 3],'sliding',@(A)any(A==labelBreast)&any(A==labelBackground)); % @(A)any(A==3)&any(A==1));
    else
        STATS = regionprops(mask==labelBreast, 'Centroid');
        XYcenter = [round(STATS(1).Centroid(2)), round(STATS(1).Centroid(1))]; % CGy, CGx
        
        % JM - Crop the mask to inside a circle that circumscribes the front part of the breast
        slope = 1.5; 
        frontBR_mask = zeros(size(mask));
        while (slope >= 1) && sum(frontBR_mask(:))==0
            tmpfrontBR_mask = mask; 
            [N,M] = size(tmpfrontBR_mask);
            for j = 1:M
                line = XYcenter(2)+floor((XYcenter(1)-j)*slope);
                if (line < 1)
                    line = 1;
                elseif (line > N)
                    line = N;
                end
                for i = 1:line % from  line:N 
                    tmpfrontBR_mask(i,j) = 0;
                end
            end
            frontBR_mask = tmpfrontBR_mask;
            slope = slope - 0.2;
        end
        %figure; imshow(frontBR_mask,[]); hold on; plot(XYcenter(2),XYcenter(1),'r*' )
        % Prepare binary mask
        frontBR_mask(frontBR_mask~=labelBreast) = -1;
        frontBR_mask(frontBR_mask==labelBreast) = 1;
        frontBR_mask(frontBR_mask~=1) = 0;
        frontBR_mask = findLCC(frontBR_mask, 6); % largest conected component
        % Boundary of binary mask
        bdImg=colfilt(frontBR_mask,[3 3],'sliding',@(A)any(A==1)&any(A==0)); 
        [yf,xf]=find(bdImg); 
        D = zeros(1, length(yf));
        for i=1:length(yf), 
            D(i) = pdist([XYcenter; yf(i) xf(i)]); 
        end; 
        radiiy = max(D);
        if(isempty(radiiy)),radiiy=450;end %hack by Michiel Kallenberg: sorry; didnt have the time to investigate te impact of this hack, but without it the code may crash

        % [radiiy, idr] = max(D); figure; imshow(frontBR_mask,[]); hold on; plot(XYcenter(2),XYcenter(1),'r*' ); plot(xf(idr), yf(idr), 'g*')
        
         % Circle
         [xx,yy]=ndgrid(1:size(mask,1), 1:size(mask,2));
         CroppingMask= uint8( (xx-XYcenter(1)).^2+(yy-XYcenter(2)).^2<=radiiy^2 );
         % MC = uint8(CroppingMask)+mask; figure; imshow(MC,[])

         mask = bsxfun(@times,  uint8(mask), uint8(CroppingMask)); 
         mask(mask==0) = labelBackground;
         %figure; imshow(mask,[])

         boundaryImg=colfilt(mask,[3 3],'sliding',@(A)any(A==labelBreast)&any(A==labelBackground)); % @(A)any(A==3)&any(A==1));        

        % exclude straight border line generated from cropped images
        idxc = find(sum(mask,1),1,'first');
        idxr = find(sum(mask,2),1,'first');
        boundaryImg(1:idxr+1,:)=0;
        boundaryImg(:,1:idxc+1)=0;   
    end
    [y,x]=find(boundaryImg);
    %MC = uint8(CroppingMask)+mask; figure; imshow(MC,[])
end

% Trying to crop Ellipse, but it needs different coord.? 
% radiix = min(pdist([ XYcenter; [linT xf(idlinT)]]), pdist([XYcenter; [ yf(idcolB) colB]]));
%         % figure; imshow(bdImg,[]); hold on; plot(XYcenter(2),XYcenter(1),'r*' )
%         % plot(xf(idlinT), linT, 'g*'); plot(colB, yf(idcolB), 'g*')
%          r_sq = [radiiy, radiix] .^ 2;  %# Ellipse radii squared (y-axis, x-axis)
%          [mX, mY] = meshgrid(1:size(mask, 2), 1:size(mask, 1));
%          ellipse_mask = (r_sq(2) * (mX - XYcenter(2)) .^ 2 + ...
%          r_sq(1) * (mY - XYcenter(1)) .^ 2 <= prod(r_sq));
%          M= uint8(ellipse_mask)+mask; figure;  imshow(M,[])


% Old running version - cropping circle based on distance to the top of pectoral line
%         yRef = min(y);
%         xRef = min(x(y==yRef));
%         Radius= pdist([XYcenter; yRef xRef]);
%         [xx,yy]=ndgrid(1:size(mask,1), 1:size(mask,2));
%         CroppingMask= uint8( (xx-XYcenter(1)).^2+(yy-XYcenter(2)).^2<=Radius^2 );
%         %mask = mask.*CroppingMask; 
%         %imshow(mask,[])
        

% exclude line from cropped images
%     [maxr, maxc] =size(mask); idxr=false(1,maxr); idxc=false(1,maxc); 
%     for r=1:maxr, idxr(r)=sum(mask(r,:)==labelBackground)==maxc; end
%     for r=1:maxc, idxc(r)=sum(mask(:,r)==labelBackground)==maxr; end
%     idxr = idxr | [idxr(2:end) 1] | [1 idxr(1:end-1)];
%     idxc = idxc | [idxc(2:end) 1] | [1 idxc(1:end-1)];
%     boundaryImg(idxr,:)=0;
%     boundaryImg(:,idxc)=0;


% Trying to cut mask based on two tangent 
% % Top
%         idx = (x>XYcenter(2));
%         idxtop = (y<XYcenter(1)*0.7);
%         xt  = x(idx & idxtop);
%         yt = y(idx & idxtop);
%         coeffsx = polyfit(xt, yt, 1);
%         fittedX = linspace(1, max(xt), 1000);
%         fittedY = polyval(coeffsx, fittedX);
%         plot(fittedX, fittedY, 'b-', 'LineWidth', 3);
%         
%         % Bottom
%         idx = (x>XYcenter(2)); % col
%         idxtop = (y>XYcenter(1)*1.5);
%         xt  = x(idx & idxtop);
%         yt = y(idx & idxtop);
%         coeffsx = polyfit(xt, yt, 1);
%         fittedX = linspace(1, max(xt), 1000);
%         fittedY = polyval(coeffsx, fittedX);
%         plot(fittedX, fittedY, 'b-', 'LineWidth', 3);
%         
%         idx = (x>XYcenter(2)); % col
%         yt = y(idx); 
%         xt  = x(idx);
%         for i=1:length(yt), D(i) = pdist([XYcenter; yt(i) xt(i)]); end; max(D)
%         figure; plot(D)
