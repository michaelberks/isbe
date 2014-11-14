%Staal Script
d_root = [asymmetryroot 'data/retinograms/dRIVE/training/'];
mkdir([d_root 'staal\orig\ridges\']);
mkdir([d_root 'staal\orig\convex_sets\']);
mkdir([d_root 'staal\orig\g2d_ori\']);
mkdir([d_root 'staal\orig\g2d_lambda\']);

for ii = 21:40
        
    %load ret and take green channel
    ret = load_uint8([d_root 'images_extended/' zerostr(ii,2) '_training_ext']);
    ret = ret(:,:,2);

    f_mask = u_load([d_root 'foveal_masks/' zerostr(ii,2) '_training_f_mask.mat']);
    
    %compute 1st deriv.
    [g1d_r g1d_o] = gaussian_1st_derivative_gradient2(ret, 1.5);

    %compute 2nd deriv. orientations at 1.5
    [g2d_r g2d_o] = gaussian_2nd_derivative_line(ret, 1.5);

    %compute second deriv responses at sigma = 0.5, 1, 2, 4
    G_r = cell(4,1);
    for jj = 1:4
        G_r{jj} = gaussian_2nd_derivative_line(ret, 2^(jj-2));
    end
    
    %Compute ridges
    ridge_map = mb_curvature_sign_change(g1d_r, g2d_r, g2d_o, 0, 0) > 0;
    ridge_map(~f_mask) = 0;
    
    %Compute convex sets
    group_map = staal_ridge_group(ridge_map, g2d_o);
    
    %Compute convex set regions
    [py px] = find(f_mask & ~group_map);
    [ry rx] = find(group_map > 0);
    set_nums = griddata(rx, ry, group_map(group_map > 0), px, py, 'nearest');    
    set_map = group_map;
    set_map(f_mask & ~group_map) = set_nums;
    
    %save outputs
    save([d_root 'staal\orig\g2d_ori\' zerostr(ii,2) '_ori_map.mat'], 'g2d_o');
    save([d_root 'staal\orig\g2d_lambda\' zerostr(ii,2) '_G_lambda.mat'], 'G_r');
    save([d_root 'staal\orig\ridges\' zerostr(ii,2) '_ridge_map.mat'], 'ridge_map');
    save([d_root 'staal\orig\convex_sets\' zerostr(ii,2) '_set_map.mat'], 'set_map');
    
    %figure; 
    %subplot(1,2,1); imgray(ridge_map);
    %subplot(1,2,2); imgray(set_map);
end
%%
%Now extract features at ridges
d_root = [asymmetryroot 'data/retinograms/dRIVE/test/'];
for ii = 1%:20
    
    
    load([d_root 'staal\orig\g2d_ori\' zerostr(ii,2) '_ori_map.mat'], 'g2d_o');
    load([d_root 'staal\orig\g2d_lambda\' zerostr(ii,2) '_G_lambda.mat'], 'G_r');
    load([d_root 'staal\orig\ridges\' zerostr(ii,2) '_ridge_map.mat'], 'ridge_map');
    load([d_root 'staal\orig\convex_sets\' zerostr(ii,2) '_set_map.mat'], 'set_map');
    
    convex_set_map = set_map;
    convex_set_map(~ridge_map) = 0;
    
    num_sets = max(convex_set_map(:));
    
    prof_width = 15;
    prof_offsets = -prof_width:prof_width;
    
    %Convex properties
    X_convex = zeros(num_sets, 18);
    set_ends = zeros(num_sets,1);
    for ss = 1:10%num_sets
        this_set_map = convex_set_map == ss;
        this_set_ends = bwmorph(this_set_map, 'endpoints');
        set_ends(ss) = sum(this_set_ends(:));
        [y x] = find(this_set_ends);
        figure; imgray(this_set_map);
        plot(x, y, 'r.');
        
    end
end
%%
        [sy sx] = find(this_set_map);
        s_theta = g2d_o(this_set_map);
        
        num_pts = length(sy);
        
        prof_x = bsxfun(@plus, prof_offsets*cos(s_theta), sx); 
        prof_y = bsxfun(@plus, -prof_offsets*sin(s_theta), sx);
        
        prof_green = interp2(ret, prof_x, prof_y, 'bilinear*');
        prof_green = mean(prof_green);
        
        %1)h
        X_convex(ss,1) = prof_green(prof_width+1);
        %2)w
        [re re_n] = min(diff(prof_green(1:prof_width+1)));
        [le le_n] = max(diff(prof_green(prof_width+1:end)));
        X_convex(ss,2) = prof_width - re_n + le_n;
        %3)h_w
        X_convex(ss,3) = h/w;
        %4)se
        X_convex(ss,4) = re + le;
        %5)se/w
        X_convex(ss,5) = se/w;
        %6)he
        X_convex(ss,6) = (prof_green(re_n) + prof_green(prof_width+le_n))/2;
        %7)h/he
        X_convex(ss,7) = h - he;
        %8)
        X_convex(ss,1) = h / he;
        
    end
        
        
end