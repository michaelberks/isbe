function new_shape = interp_shape(shape_pts, method)

    if nargin < 2
        method = 'linear';
    end

    [dummy unique_idx] = unique(shape_pts, 'rows');
    shape_pts = shape_pts(sort(unique_idx), :);

    s_diff = diff(shape_pts);
    s_lengths = [0; cumsum(sqrt(s_diff(:,1).^2 + s_diff(:,2).^2))];
    clear s_diff;
    
    new_shape = interp1(s_lengths, shape_pts, 0:s_lengths(end), method);