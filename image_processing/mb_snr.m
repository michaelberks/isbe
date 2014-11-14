function [snr local_sd local_mean] = mb_snr(im_in, win_size)

    im_in = double(im_in);
    local_mean = imfilter(im_in, ones(win_size)/(win_size^2), 'symmetric');
    local_sq = imfilter(im_in.^2, ones(win_size)/(win_size^2), 'symmetric');
    local_sd = sqrt(local_sq - (local_mean.*local_mean));
   
    
    snr = local_mean ./ local_sd;
end