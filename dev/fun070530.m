function increase = fun070530(masses, n1, n2, spacing, sigma)
    
    increase(length(masses)) = 0;
    for ii = 1:length(masses);
        temp = load(['annotations\', masses(ii).name]);
        mass = temp.mass; clear temp;
        
        [r c] = size(mass.mass_ROI);
        mass_bw = poly2mask(mass.mass_outline(:,1),...
                            mass.mass_outline(:,2),...
                            r, c);

        %dilate the shape mask by n1
        mask1 = mass_bw;
        for jj = 1:n1
            mask1 = imdilate(mask1, strel('disk', 1));
        end

        %dilate the shape mask by n2
        mask2 = mask1;
        for jj = 1:n2
            mask2 = imdilate(mask2, strel('disk', 1));
        end

        %subtract mask1 from mask2, result is ring in shape of mass, n1 pixels
        % from the shape border, n2-n1 pixels wide
        mask3 = mask2 - mask1; clear mask2 mask1;

        %create a mask of equally spaced dots (as defined by input: spacing) 
        mask4 = zeros([r c]);
        mask4(1:spacing:end, 1:spacing:end) = 1;

        %mask of landmark pts is intersection of spaced dots and ring of shape
        mask5 = mask3 & mask4; clear mask3 mask4
        landmark_pts = find(mask5); clear mask5;

        smooth_ROI = wiener2(mass.mass_ROI, [sigma sigma]);
        increase(ii) = (mean(smooth_ROI(landmark_pts)) - ...
            mean(mass.mass_ROI(landmark_pts)));
        clear smooth_ROI;
    end
end