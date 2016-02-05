function nipple = locateNippleParallelPectoral(mask, lineprm, labelBreast, labelBackground)

    % Restrict the nipple to the area in front of the center of gravity 
    restrictNippleJM = 1;
    restrinctMask = 1;
    
    if restrictNippleJM || restrinctMask
        % the center of gravity (CG)
        STATS = regionprops(mask==labelBreast, 'Centroid');
        CG = round(STATS(1).Centroid); % CGy, CGx

        if restrinctMask
            [ ~, ~, mask] = extractCroppedBoundary(mask,labelBreast, labelBackground); % , [CG(2) CG(1)] xp(2),xp(1));
        end
    end
    
    [N,M] = size(mask);
    nipple = [];
    for line_rho = lineprm(1):sqrt(N^2+M^2)
        for I = N:-1:1
            J = round((line_rho - (I-0.5) * sin(lineprm(2))) / ...
                cos(lineprm(2)) + 0.5);
            % original Sami: (J < 1)   -    Alternatively JM: (J < CG(1)) || (J > CG(1)+d)
            if ( restrictNippleJM &&  J < CG(1) ) ...  % || (J > CG(1)+d))
                    || ( ~restrictNippleJM && (J < 1) ) 
                continue;
            elseif J > M
                break;
            elseif mask(I,J) == labelBreast
                nipple = [I,J];
                break;
            end            
        end
    end
end


%         if restrict Nipple by distance of Centre gravity
%             % Intersection  with pectoral line perpendicular and the center of gravity (CG)
%             im_crop = [0 0 0 0 ];
%             [N,M] = size(mask);
%             pectoral = LineConv(convertPolarPectLine(mask, lineprm, im_crop)); % x, y
%             xp=cross(cross([pectoral(1);pectoral(2);0],[CG(1); CG(2) ;1]), pectoral); 
%             xp=reshape(xp/xp(3),3,1);
%             d = pdist([CG(1) CG(2); xp(1) xp(2)],'euclidean'); % x, y        
%         end