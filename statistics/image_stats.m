function [local_sd im_mean local_sd_map] = image_stats(im_in, win_size)

    im_in = double(im_in);
    t1 = imfilter(im_in, ones(win_size)/(win_size^2), 'symmetric');
    t2 = imfilter(im_in.^2, ones(win_size)/(win_size^2), 'symmetric');
    local_sd_map = t2 - (t1.*t1);
   
    local_sd = mean(local_sd_map(:));
    im_mean = mean(im_in(:));
    if nargout > 2
        local_sd_map = local_sd_map - min(local_sd_map(:));
        local_sd_map = local_sd_map / max(local_sd_map(:));
    end
end