function lp = fromPolar2LineConv(N, M, lineprm, im_crop)

    [x1,y1,x2,y2] = polar2borderpoints([N M],lineprm(1),lineprm(2));
    
    % Adjust points
    if numel(lineprm) > 0
        x1 = x1 + im_crop(1);
        x2 = x2 + im_crop(1);

        y1 = y1 + im_crop(3);
        y2 = y2 + im_crop(3);

        lp = [x1,y1;x2,y2];
    else
        lp = [];
    end
    lp = LineConv(lp);
    lp = lp/lp(3);
end