
orientations = [Inf -Inf; 0 pi; pi 2*pi; 2*pi 3*pi; -3*pi -2*pi; -2*pi -pi; -pi 0] / 6;
for ii = 1:3

    %load tree
    dt = u_load(['C:\isbe\dev\background\dual_tree\mass_2\mass', zerostr(ii,3) ,'_dual_tree.mat']);
    
    for jj = 1:5; dt{jj} = real(dt{jj}); end
    
    %Load CLS result for each mass
    cls = u_load(['C:\isbe\dev\background\cls\mass_2\mass_cls_', zerostr(ii,3), '.mat']);
    
    %Generate oriented CLS map
    [cls_map] =  mb_cls_map_orientations(cls, orientations);
    cls_lev = imresize(cls_map{1}, size(dt{1}(:,:,1)), 'nearest');
    
    figure; hold on;

    [rows cols] = size(dt{1}(:,:,1));

    for r = 1:rows
        for c = 1:cols;
            
            if ~cls_lev(r, c)
                %Non-CLS
                plot3(dt{1}(r,c,1), dt{2}(ceil(r/2),ceil(c/2),1), ...
                    dt{3}(ceil(r/4),ceil(c/4),1), 'rx');
            elseif cls_lev(r, c) == 2
                %CLS aligned
                plot3(dt{1}(r,c,1), dt{2}(ceil(r/2),ceil(c/2),1), ...
                    dt{3}(ceil(r/4),ceil(c/4),1), 'gx');
            else
                %CLS but non-aligned
                plot3(dt{1}(r,c,1), dt{2}(ceil(r/2),ceil(c/2),1), ...
                    dt{3}(ceil(r/4),ceil(c/4),1), 'bx');
            end
        end
    end
end

%%
orientations = [Inf -Inf; 0 pi; pi 2*pi; 2*pi 3*pi; -3*pi -2*pi; -2*pi -pi; -pi 0] / 6;

for ii = 1:10 %53 = number of masses
    

    %Load CLS result for each mass
    cls = u_load(['C:\isbe\dev\background\cls\mass_2\mass_cls_', zerostr(ii,3), '.mat']);
    
    
    %Generate oriented CLS map
    [cls_map] =  mb_cls_map_orientations(cls, orientations);
    
    %Load pyramid
    dt = u_load(['C:\isbe\dev\background\dual_tree\mass_2\mass', zerostr(ii,3) ,'_dual_tree.mat']);
    
    for jj = 1:5; dt{jj} = real(dt{jj}); end
    
    cls_1 = imresize(cls_map{1}, size(dt{1}(:,:,1)), 'nearest');
    
    [r1c c1c] = find(cls_1);
    [r1n c1n] = find(~cls_1);
    
    idx1c = sub2ind(size(dt{1}(:,:,1)), r1c, c1c);
    idx1n = sub2ind(size(dt{1}(:,:,1)), r1n, c1n);
    
    r2c = ceil(r1c/2); c2c = ceil(c1c/2);
    r2n = ceil(r1n/2); c2n = ceil(c1n/2);
    
    idx2c = sub2ind(size(dt{2}(:,:,1)), r2c, c2c);
    idx2n = sub2ind(size(dt{2}(:,:,1)), r2n, c2n);
    
    r3c = ceil(r1c/4); c3c = ceil(c1c/4);
    r3n = ceil(r1n/4); c3n = ceil(c1n/4);
    
    idx3c = sub2ind(size(dt{3}(:,:,1)), r3c, c3c);
    idx3n = sub2ind(size(dt{3}(:,:,1)), r3n, c3n);  
    
    figure;
    
    for ori = 1:6           

        temp1 = reshape(dt{1}(:,:,ori),[],1);
        temp2 = reshape(dt{2}(:,:,ori),[],1);
        temp3 = reshape(dt{3}(:,:,ori),[],1);
        
        data_non_cls = [temp1(idx1n) temp2(idx2n) temp3(idx3n)];
                    
        data_cls = [temp1(idx1c) temp2(idx2c) temp3(idx3c)];
        
        [P, L, B, mu] = st_pca(data_non_cls, 1);
        
        [x y z] = mb_ellipsoid(mu, 1.96*sqrt(L), P);
         
        subplot(2,3,ori); hold on;
        
        mesh(x, y, z, 'FaceColor', 'none');
        
        [P, L, B, mu] = st_pca(data_cls, 1);
        [x y z] = mb_ellipsoid(mu, 1.96*sqrt(L), P);
        
        mesh(x, y, z, 'FaceColor', 'none');
        
    end
    
    clear cls cls_map cls_lev dt pyr temp* x y z
end
%%
for ii = 1:10 %53 = number of masses
    

    %Load CLS result for each mass
    cls = u_load(['C:\isbe\dev\background\cls\mass_2\mass_cls_', zerostr(ii,3), '.mat']);
    
    
    %Generate oriented CLS map
    [cls_map] =  mb_cls_map_orientations(cls);
    
    %Load pyramid
    pyr = u_load(['C:\isbe\dev\background\pyramid\mass_2\mass', zerostr(ii,3) ,'_pyramid.mat']);
    
%     [r1c c1c] = find(cls_map{1});
%     [r1n c1n] = find(~cls_map{1});
%     
%     idx1c = sub2ind(size(pyr{2,1}), r1c, c1c);
%     idx1n = sub2ind(size(pyr{2,1}), r1n, c1n);
%     
%     r2c = ceil(r1c/2); c2c = ceil(c1c/2);
%     r2n = ceil(r1n/2); c2n = ceil(c1n/2);
%     
%     idx2c = sub2ind(size(pyr{3,1}), r2c, c2c);
%     idx2n = sub2ind(size(pyr{3,1}), r2n, c2n);
%     
%     r3c = ceil(r1c/4); c3c = ceil(c1c/4);
%     r3n = ceil(r1n/4); c3n = ceil(c1n/4);
%     
%     idx3c = sub2ind(size(pyr{4,1}), r3c, c3c);
%     idx3n = sub2ind(size(pyr{4,1}), r3n, c3n);  
    
%     figure;
    
    display(['mass ', num2str(ii)]);
    for ori = 1:4
        
        
        [r1c c1c] = find(cls_map{1} == ori);
        [r1n c1n] = find(~cls_map{1});
        
        data_cls = zeros(size(r1c,1), 5);
        data_non_cls = zeros(size(r1n,1), 5);

        idx_c = sub2ind(size(pyr{2,1}), r1c, c1c);
        idx_n = sub2ind(size(pyr{2,1}), r1n, c1n);
        
        data_cls(:,1) = pyr{2,ori}(idx_c);
        data_non_cls(:,1) = pyr{2,ori}(idx_n);
        
        for level = 2:5
            rc = ceil(r1c/2^(level-1));
            cc = ceil(c1c/2^(level-1));
            rn = ceil(r1n/2^(level-1));
            cn = ceil(c1n/2^(level-1));

            idx_c = sub2ind(size(pyr{level+1,1}), rc, cc);
            idx_n = sub2ind(size(pyr{level+1,1}), rn, cn);
            
            data_cls(:,level) = pyr{level+1,ori}(idx_c);
            data_non_cls(:,level) = pyr{level+1,ori}(idx_n);

 

        end
        
        [P, L, B, mu] = st_pca(data_non_cls, 0.95);
        display(['Non CLS modes: ', num2str(size(P,2))]);
        
%         [x y z] = mb_ellipsoid(mu, 1.96*sqrt(L), P);
%          
%         subplot(2,2,ori); hold on;
%         
%         mesh(x, y, z, 'FaceColor', 'none');
        
        [P, L, B, mu] = st_pca(data_cls, 0.95);
        display(['CLS modes: ', num2str(size(P,2))]);
        display('');
%         [x y z] = mb_ellipsoid(mu, 1.96*sqrt(L), P);
%         
%         mesh(x, y, z, 'FaceColor', 'none');
        
    end
    
    clear cls cls_map cls_lev dt pyr temp* x y z
end

%%
for ii = 1:10 %53 = number of masses
    

    %Load CLS result for each mass
    cls = u_load(['C:\isbe\dev\background\cls\mass_2\mass_cls_', zerostr(ii,3), '.mat']);
    
    
    %Generate oriented CLS map
    [cls_map] =  mb_cls_map_orientations(cls);
    
    %Load pyramid
    pyr = u_load(['C:\isbe\dev\background\pyramid\mass_2\mass', zerostr(ii,3) ,'_pyramid.mat']);
    
    [r1c c1c] = find(cls_map{1} > 0);
    [r1n c1n] = find(~cls_map{1});

    data_cls = zeros(size(r1c,1), 20);
    data_non_cls = zeros(size(r1n,1), 20);
    
    display(['mass ', num2str(ii)]);
    
    for level = 1:5 
        
        rc = ceil(r1c/2^(level-1));
        cc = ceil(c1c/2^(level-1));
        rn = ceil(r1n/2^(level-1));
        cn = ceil(c1n/2^(level-1));

        idx_c = sub2ind(size(pyr{level+1,1}), rc, cc);
        idx_n = sub2ind(size(pyr{level+1,1}), rn, cn);
        
        for ori = 1:4
            
            idx = (level - 1)*4 + ori;
            data_cls(:,idx) = pyr{level+1,ori}(idx_c);
            data_non_cls(:,idx) = pyr{level+1,ori}(idx_n);
        end 
    end
    
    [P, L, B, mu] = st_pca(data_non_cls, 0.95);
    display(['Non CLS modes: ', num2str(size(P,2))]);

    [P, L, B, mu] = st_pca(data_cls, 0.95);
    display(['CLS modes: ', num2str(size(P,2))]);
    display([]);
    
    clear cls cls_map cls_lev dt pyr temp* x y z
end

%%
for ii = 1:10 %53 = number of masses
    

    %Load CLS result for each mass
    cls = u_load(['C:\isbe\dev\background\cls\mass_2\mass_cls_', zerostr(ii,3), '.mat']);
    
    
    %Generate oriented CLS map
    [cls_map] =  mb_cls_map_orientations(cls);
    
    %Load pyramid
    dt = u_load(['C:\isbe\dev\background\dual_tree\mass_2\mass', zerostr(ii,3) ,'_dual_tree.mat']);
    
    [r1c c1c] = find(cls_map{1} > 0);
    [r1n c1n] = find(~cls_map{1});

    data_cls = zeros(size(r1c,1), 60);
    data_non_cls = zeros(size(r1n,1), 60);
    
    display(['mass ', num2str(ii)]);
    
    for level = 1:5 
        
        rc = ceil(r1c/2^(level));
        cc = ceil(c1c/2^(level));
        rn = ceil(r1n/2^(level));
        cn = ceil(c1n/2^(level));

        idx_c = sub2ind(size(dt{level}(:,:,1)), rc, cc);
        idx_n = sub2ind(size(dt{level}(:,:,1)), rn, cn);
        
        for ori = 1:6
            
            idx_real = (level - 1)*4 + ori;
            idx_imag = (level - 1)*4 + ori + 30;
            
            temp = real(dt{level}(:,:,ori));
            data_cls(:,idx_real) = temp(idx_c);
            data_non_cls(:,idx_real) = temp(idx_n);
            
            temp = imag(dt{level}(:,:,ori));
            data_cls(:,idx_imag) = temp(idx_c);
            data_non_cls(:,idx_imag) = temp(idx_n);
            clear temp
        end 
    end
    
    [P, L, B, mu] = st_pca(data_non_cls, 0.95);
    display(['Non CLS modes: ', num2str(size(P,2))]);

    [P, L, B, mu] = st_pca(data_cls, 0.95);
    display(['CLS modes: ', num2str(size(P,2))]);
    display([]);
    
    clear cls cls_map cls_lev dt pyr temp* x y z
end