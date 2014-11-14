function error_struct = detection_errors(imlists)

if (nargin==0)
    % dummy mode - return empty error structure
    error_struct = struct('Az', NaN);
    return
end

values = [];
labels = [];
for i = 1:length(imlists)
    images = load_vessel_images(imlists(i));
    values = [values; images.probability(images.foveal_mask)];
    labels = [labels; images.vessel_mask(images.foveal_mask)];
end

[roc, Az] = calculate_roc_curve(values, labels);

error_struct = struct('Az', Az);

