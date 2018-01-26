function [local_mean local_std local_x2] = image_stats(im_in, win_size)

    im_in = double(im_in);
    %t1 = imfilter(im_in, ones(win_size)/(win_size^2), 'symmetric');
    %t2 = imfilter(im_in.^2, ones(win_size)/(win_size^2), 'symmetric');
    local_x = conv2(im_in, ones(win_size)/(win_size^2), 'valid');
    local_x2 = conv2(im_in.^2, ones(win_size)/(win_size^2), 'valid');
    local_var = local_x2 - (local_x.*local_x);
   
    local_mean = local_x;
    local_std = sqrt(local_var);
    
end