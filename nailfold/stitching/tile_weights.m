function weights = tile_weights(image_size, weight_shape)

rows = image_size(1);
cols = image_size(2);

cx = cols/2;
cy = rows/2;

switch weight_shape
    case 'rect'
        %Create an array of tile weights that fades the contribution of each tile
        %horizontally from centre to edge
        weights = linspace(0, 1, ceil(cx)+1);
        weights = weights(2:end); % remove the zero weight
        if rem(cols,2)
            weights = repmat([weights weights(end-1:-1:1)], [rows, 1]);
        else
            weights = repmat([weights weights(end:-1:1)], [rows, 1]);
        end
        
    case 'circle'
        %Create an array of tile weights that fades the contribution of each tile
        %circularly from centre to edge
        x_tile = repmat(1:cols, rows, 1);
        y_tile = repmat((1:rows)', 1, cols);
        weights = sqrt( (0.5*(x_tile - cx)/cols).^2 + (0.5*(y_tile - cy)/rows).^2 );
        weights = weights / (max(weights(:)) + 1e-6);
        weights = 1 - weights;
        
    case 'diamond'
        %Create an array of tile weights that fades the contribution of each tile
        %horizontally from centre to edge
        weights = linspace(0, 1, ceil(cx)+1);
        weights = weights(2:end); % remove the zero weight
        if rem(cols,2)
            weights = repmat([weights weights(end-1:-1:1)], [rows, 1]);
        else
            weights = repmat([weights weights(end:-1:1)], [rows, 1]);
        end
        
        weights2 = linspace(0, 1, ceil(cy)+1);
        weights2 = weights2(2:end); % remove the zero weight
        if rem(rows,2)
            weights2 = repmat([weights2 weights2(end-1:-1:1)], [cols, 1])';
        else
            weights2 = repmat([weights2 weights2(end:-1:1)], [cols, 1])';
        end
        weights = min(weights, weights2);
        
    case 'uniform'
        %Create uniform weights
        weights = ones(rows, cols);
end

