function [pixel_list, this_SSD] = find_matches(sample_window, I_sample, weighting_matrix, ErrThreshold, row, col)
    
    % Taken from the Efros Pseudocode
    win_size = (size(sample_window,1) - 1);
    weighting_matrix(isnan(sample_window)) = 0;
    weighting_matrix = weighting_matrix(:) ./ sum(weighting_matrix(:));
    sample_window(isnan(sample_window)) = 0;
    
    this_SSD = repmat(Inf, size(I_sample));

    for i = 1: size(I_sample,1) - win_size
        for j = 1:size(I_sample,2) - win_size
            dist = (sample_window - I_sample(i:i+win_size, j:j+win_size)).^2;
            dist = dist(:) .* weighting_matrix;   
            this_SSD(i+(win_size/2),j+(win_size/2)) = sum(dist);
        end
    end

    this_SSD(row, col) = Inf;
    pixel_list = find(this_SSD <= min(this_SSD(:))*(1+ErrThreshold));
end