function [pixel_list, this_SSD] = find_matches2(sample_window, I_sample, weighting_matrix, ErrThreshold)
    
    % MB
    nan_idx = isnan(sample_window);
    weighting_matrix(nan_idx) = 0;
    this_SSD = mb_normxcorr2(sample_window, I_sample, weighting_matrix);
    this_SSD = padarray(this_SSD, ones(1,2)*(size(sample_window, 1)-1)/2, -2); %use -2 cos NaN buggers up sort
%    pixel_list = find(this_SSD >= max(this_SSD(:))*(1-ErrThreshold));
%     [dummy pixel_list] = max(this_SSD(:));
    [dummy pixel_list] = sort(this_SSD(:), 'descend');
end