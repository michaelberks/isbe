function [enhance] = fast_linop(im_orig, dim)

%set up original line
    rad = linspace(0, pi, 13);
    rad(end) = [];
    len = length(rad);
    [line_k1 sq_k1] = create_kernels;
    [lin_stren1 lin_ori1] = filter_response(line_k1, sq_k1);
    
    [line_k2 sq_k2] = scale_kernels(line_k1, sq_k1);
    [lin_stren2 lin_ori2] = filter_response(line_k2, sq_k2);
    
    [line_k3 sq_k3] = scale_kernels(line_k2, sq_k2);
    [lin_stren3 lin_ori3] = filter_response(line_k3, sq_k3);
    
    [line_k4 sq_k4] = scale_kernels(line_k3, sq_k3);
    [lin_stren4 lin_ori4] = filter_response(line_k4, sq_k4);
    
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
    max_stren = reshape(max_stren, size(im_orig));
    max_ori = reshape(all_ori(sub2ind(size(all_ori), ...
        scales, 1:size(im_orig(:)))), size(im_orig)); clear all_ori;
    scales = reshape(scales, size(im_orig));
    
    figure; imagesc(max_stren); axis image; colormap gray; hold on;
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
    
    function [l_k s_k] = create_kernels
        
        [x y] = meshgrid([1-dim:dim-1]);

        line_im = zeros(2*dim-1, 2*dim-1); 
        line_im(:, dim) = 1;

        sq_v = [1-dim dim-1 dim-1 1-dim; 1-dim 1-dim dim-1 dim-1];
        for ii = 1:len
            rot_mat = [cos(rad(ii)) -sin(rad(ii)); sin(rad(ii)) cos(rad(ii))];

            line_grid = rot_mat*[y(:)'; x(:)'] + dim;


            sq_grid = rot_mat*sq_v + 2*dim;

            l_k(ii).conv = zeros(2*dim-1, 2*dim-1);
            l_temp = interp2(line_im, line_grid(2,:), line_grid(1,:));

            l_temp(isnan(l_temp)) = 0;

            l_k(ii).conv(sub2ind([2*dim-1, 2*dim-1],...
                y(:)+dim, x(:)+dim)) = l_temp;
            l_k(ii).conv = l_k(ii).conv / sum(l_k(ii).conv(:));

            s_k(ii).conv = poly2mask(sq_grid(1,:), sq_grid(2,:), 4*dim-1, 4*dim-1);
            s_k(ii).conv = s_k(ii).conv / sum(s_k(ii).conv(:));
        end
    end

    function [l_stren l_ori] = filter_response(l_k, s_k)
        for ii = 1:12
            lr_ii = imfilter(im_orig, l_k(ii).conv, 'replicate');
            lr(ii,:) = lr_ii(:)';

            sqr_ii = imfilter(im_orig, s_k(ii).conv, 'replicate');
            sqr(ii,:) = sqr_ii(:)';
        end
    
        [max_lr l_ori] = max(lr);
        max_lr = reshape(max_lr, size(im_orig));
        ls = reshape(sqr(sub2ind(size(sqr), ...
            l_ori, 1:length(im_orig(:)))), size(im_orig));

        l_ori = reshape(l_ori, size(im_orig));
        l_stren = max_lr - ls;
        l_stren = normalise_im(l_stren);
        figure; imagesc(l_stren); colormap gray; axis image;
    end
    
    function [new_lk new_sk] = scale_kernels(old_lk, old_sk)
        [new_lr new_lc] = size(old_lk(1).conv);
        new_lr = 2*new_lr - 1;
        new_lc = 2*new_lc - 1;
        
        [new_sr new_sc] = size(old_sk(1).conv);
        new_sr = 2*new_sr - 1;
        new_sc = 2*new_sc - 1;
        
        for ii = 1:12
            new_lk(ii).conv = zeros(new_lr, new_lc);
            
            new_lk(ii).conv(1:2:end, 1:2:end) = old_lk(ii).conv;
            
            new_sk(ii).conv = zeros(new_sr, new_sc);
            new_sk(ii).conv(1:2:end, 1:2:end) = old_sk(ii).conv;
        end
    end

    function im_norm = normalise_im(im_un)
        im_un = double(im_un);
        im_norm = im_un - min(im_un(:)); clear im_un;
        im_norm = im_norm / max(im_norm(:));
    end
end
    
    