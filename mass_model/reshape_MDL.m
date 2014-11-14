function shapes_out = reshape_MDL(shapes_in)

if ndims(shapes_in) == 2
    shapes_out = permute(reshape(shapes_in',...
        [], 2, size(shapes_in, 1)), [2 1 3]);
elseif ndims(shapes_in) == 3
    shapes_out = reshape(permute(shapes_in,...
        [2 1 3]), [], size(shapes_in, 3))';
else
    error('Shape matrix must have dimension 2 or 3');
end
    
    