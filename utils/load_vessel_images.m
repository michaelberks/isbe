function images = load_vessel_images(imlist)
% Create a structure that holds the filenames for all images of interest to
% a vessel predictor

% Since the images are in all different formats, we'll get this warning a
% lot so disable it temporarily
warning('off','load_uint8:missing_variables');

names = fieldnames(imlist);
for name = names(:)'
    fieldname = name{1};
    if ~isempty(imlist.(fieldname))
        images.(fieldname) = load_uint8(imlist.(fieldname));
    end
end

% Reenable warnings
warning('on','load_uint8:missing_variables');

