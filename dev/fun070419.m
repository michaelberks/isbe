function fun070419(roi, lev)

    cm = contourc(double(roi), [lev lev]);
    pts_vec = [];
    idx_vec = [];
    idx = 1;
    while idx < size(cm, 2)
        idx_vec(end+1) = idx;
        pts_vec(end+1) = cm(2, idx);
        idx = idx + cm(2, idx) + 1;
    end
    
    [m idx] = max(pts_vec);
    idx = idx_vec(idx);
    x_contour = cm(1, idx+1:5:idx+m);
    y_contour = cm(2, idx+1:5:idx+m);
    
    figure; imagesc(roi); axis image; hold on;
    plot(x_contour, y_contour, 'k');
end
