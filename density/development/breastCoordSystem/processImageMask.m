function imageMask = processImageMask(imageMask)

    if (isstruct(imageMask) && isfield(imageMask,'Format') && strcmp(imageMask.Format,'DICOM'))
        imageMask = dicomread(imageMask);
        imageMask = 1 - log(single(imageMask));
    end        
        
end