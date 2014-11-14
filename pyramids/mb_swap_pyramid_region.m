function new_pyr = mb_swap_pyramid_region(p1, p2, roi, subbands, levels, orientations)

%p1 = target pyramid
%p2 = source pyramid

%roi = indices of region to swap in source pyramid

%subbands = subbands at which to make the swap (if blank do all). Specified
%a nx2 array where each row gives the level and orientation of a specific
%subband. Leave empty to use level/orientation range

%levels = Vector of levels to swap at. Used to populate subbands automatically

%orientations = Vectors of orientations to swap at. Used to populate
%subbands automatically

%check if subbands needs to populated automatically
if isempty(subbands)
    [levs oris] = meshgrid(levels(1):levels(2), orientations(1):orientations(2));
    subbands = [levs(:) oris(:)];
end

new_pyr = p1; clear p1;
[rows cols] = ind2sub(size(p2{2, 1}), roi);

%Swap region at each subband
for ii = 1:size(subbands, 1)
    
    level = subbands(ii,1);
    ori = subbands(ii, 2);
    
    %check if we need to downsample indices (note we could cache the
    %downsampled indices so we don't recalculate for each orientation. If
    %the function is to be used regularly (i.e. as part of a function in
    %big for loop) it will probably be a good idea to do this)
    
    if level > 2
        new_rows = ceil(rows/(2^(level-2))); 
        new_cols = ceil(cols/(2^(level-2)));

        new_roi1 = unique(sub2ind(size(new_pyr{level, 1}), new_rows, new_cols));
        new_roi2 = unique(sub2ind(size(p2{level, 1}), new_rows, new_cols));
    else
        new_roi1 = sub2ind(size(new_pyr{level, 1}), rows, cols);
        new_roi2 = roi;
    end
    
    new_pyr{level, ori}(new_roi1) = p2{level, ori}(new_roi2);
    
end
    
    




