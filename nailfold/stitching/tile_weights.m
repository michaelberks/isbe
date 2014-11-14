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
        weights = (x_tile - cx).^2 + (y_tile - cy).^2;
        weights = weights / (max(weights(:)) + 1e-6);
        weights = 1 - weights;
        
    case 'uniform'
        %Create uniform weights
        weights = ones(rows, cols);
end

