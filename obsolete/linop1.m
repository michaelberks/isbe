function [enhance] = linop1(orig_im, dim)

    figure; imagesc(orig_im); axis image; colormap gray; hold on;
    [rows cols] = size(orig_im);
    
    rad = linspace(0, pi, 12);
    len = length(rad);
    c = cos(rad);
    s = sin(rad);
    rot = [reshape([c; -s], 1, []); reshape([s; c], 1, [])]'; clear c s;
        
    gf = zeros(2*dim + 1, 2*dim + 1);
    gf(dim+1, dim+1) = 1;        
    im_scale1 = imfilter(orig_im, gf, 'replicate', 'full');
    im_in = im_scale1;
    scale = 1;
    v_ori = (scale*(dim-1)*rot*[0 0; -1, 1]);
    [vx vy] = meshgrid(scale*[1-dim:dim-1]);
    v_stren(1,:) = vx(:); clear vx;
    v_stren(2,:) = vy(:); clear vy;
    [X Y] = meshgrid(dim+1:rows+dim, dim+1:cols+dim);
    [lin_ori1 lin_stren1] = arrayfun(@filter_response, Y, X);
    lin_stren1 = normalise_im(lin_stren1);
    clear v_stren;
    
    gf = fspecial('gaussian', 2*dim + 1, 1);
    im_scale2 = imfilter(im_scale1, gf, 'replicate', 'full');
    im_in = im_scale2;
    scale = 2;
    v_ori = (scale*(dim-1)*rot*[0 0; -1, 1]);
    [vx vy] = meshgrid(scale*[1-dim:dim-1]);
    v_stren(1,:) = vx(:); clear vx;
    v_stren(2,:) = vy(:); clear vy;
    [X Y] = meshgrid(2*dim+1:rows+2*dim, 2*dim+1:cols+2*dim);
    [lin_ori2 lin_stren2] = arrayfun(@filter_response, Y, X);
    lin_stren2 = normalise_im(lin_stren2);
    clear v_stren;
    
    im_scale3 = imfilter(im_scale2, gf, 'replicate', 'full');
    im_in = im_scale3;
    scale = 3;
    v_ori = (scale*(dim-1)*rot*[0 0; -1, 1]);
    [vx vy] = meshgrid(scale*[1-dim:dim-1]);
    v_stren(1,:) = vx(:); clear vx;
    v_stren(2,:) = vy(:); clear vy;
    [X Y] = meshgrid(3*dim+1:rows+3*dim, 3*dim+1:cols+3*dim);
    [lin_ori3 lin_stren3] = arrayfun(@filter_response, Y, X);
    lin_stren3 = normalise_im(lin_stren3);
    clear v_stren;
    
    im_scale4 = imfilter(im_scale3, gf, 'replicate', 'full');
    im_in = im_scale4;
    scale = 4;
    v_ori = (scale*(dim-1)*rot*[0 0; -1, 1]);
    [vx vy] = meshgrid(scale*[1-dim:dim-1]);
    v_stren(1,:) = vx(:); clear vx;
    v_stren(2,:) = vy(:); clear vy;
    [X Y] = meshgrid(4*dim+1:rows+4*dim, 4*dim+1:cols+4*dim);
    [lin_ori4 lin_stren4] = arrayfun(@filter_response, Y, X);
    lin_stren4 = normalise_im(lin_stren4);
    clear v_stren;
    
    all_scale(1,:) = lin_stren1(:)';
    all_scale(2,:) = lin_stren2(:)';
    all_scale(3,:) = lin_stren3(:)';
    all_scale(4,:) = lin_stren4(:)';
    
    sum_stren = lin_stren1 + lin_stren2 + lin_stren3 + lin_stren4;
    sum_stren = normalise_im(sum_stren);
    
    all_ori(1,:) = lin_ori1(:)';
    all_ori(2,:) = lin_ori2(:)';
    all_ori(3,:) = lin_ori3(:)';
    all_ori(4,:) = lin_ori4(:)';
    
    [max_stren scales] = max(all_scale); clear all_scale;
    max_stren = reshape(max_stren, rows, cols);
    max_ori = reshape(all_ori(sub2ind(size(all_ori), ...
        scales, 1:rows*cols)), rows, cols); clear all_ori;
    scales = reshape(scales, rows, cols);
    
    %figure; imagesc(lin_stren1); axis image; colormap gray; hold on;
    %figure; imagesc(lin_stren2); axis image; colormap gray; hold on;
    %figure; imagesc(lin_stren3); axis image; colormap gray; hold on;
    %figure; imagesc(lin_stren4); axis image; colormap gray; hold on;
    %figure; imagesc(max_stren); axis image; colormap gray; hold on;
    figure; imagesc(sum_stren); axis image; colormap gray; hold on;
    
    enhance.strength1 = lin_stren1; clear lin_stren1;
    enhance.strength2 = lin_stren2; clear lin_stren2;
    enhance.strength3 = lin_stren3; clear lin_stren3;
    enhance.strength4 = lin_stren4; clear lin_stren4;
    
    enhance.ori1 = lin_ori1; clear lin_ori1;
    enhance.ori2 = lin_ori2; clear lin_ori2;
    enhance.ori3 = lin_ori3; clear lin_ori3;
    enhance.ori4 = lin_ori4; clear lin_ori4;
    
    enhance.sum_strengths = sum_stren; clear sum_stren;
    enhance.max_strengths = max_stren; clear max_stren;
    enhance.max_ori = max_ori; clear max_ori;
    enhance.scales = scales; clear scales;
    

    function [ori lin_stren] = filter_response(row, col)
        [ori max_r] = line_orientation(col, row);
        lin_stren = line_strength(col, row, ori, max_r);
    end
    
    function [ls] = line_strength(x, y, rad, max_r)
        rot_mat = [cos(rad) -sin(rad); sin(rad) cos(rad)];
        v = rot_mat*v_stren;
        v(1,:) = v(1,:) + x;
        v(2,:) = v(2,:) + y;
        ls = max_r - ...
            sum(interp2(double(im_in), v(1,:), v(2,:), '*linear'))/(2*dim -1)^2;
        
    end

    function [ori max_r] = line_orientation(x, y)
        
        v = v_ori + repmat([x x; y y], len, 1);
        for ii = 1:len
            jj = 2*ii;
            response(ii) = sum(improfile(im_in, v(jj-1,:), v(jj,:), 2*dim-1));
        end
        [max_r idx_r] = max(response);
        max_r = max_r/(2*dim-1);
        ori = rad(idx_r);
    end

    function im_norm = normalise_im(im_un)
        im_norm = im_un - min(im_un(:)); clear im_un;
        im_norm = im_norm / max(im_norm(:));
    end
end