function pectoral = convertPolarPectLine(mask, lineprm, im_crop)
    [N,M] = size(mask);
    if numel(lineprm) > 0
        [x1,y1,x2,y2] = polar2borderpoints([N M],lineprm(1),lineprm(2));
    end

    % Adjust points
    if numel(lineprm) > 0
        x1 = x1 + im_crop(1);
        x2 = x2 + im_crop(1);
        y1 = y1 + im_crop(3);
        y2 = y2 + im_crop(3);
        pectoral = [x1,y1;x2,y2];
    else
        pectoral = [];
    end

end