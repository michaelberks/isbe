load C:\isbe\matlab_code\trunk\wavelet_dtcwt_toolbox4_1\lenna
dt = dtwavexfm2(X, 5);
[ilp icp] = mb_dual_tree_transform(dt);
%
%figure;
for level = 1:4
    dim = 2^level;
    st = (dim+1)/2;
    [max_icp] = max(icp{level}, [], 3);

    figure;
    imagesc(X); axis image; colormap(gray(256)); hold on;
    quiver(st:dim:257-st, st:dim:257-st, real(max_icp), -imag(max_icp), 2, 'r'); 
    title(['Level ', num2str(level), ' ICP coefficients']); 
    
end
%%
for level = 1:4
    dim = 2^level;
    st = (dim+1)/2;
    [max_ilp] = max(ilp{level}, [], 3);

    figure;
    imagesc(X); axis image; colormap(gray(256)); hold on;
    quiver(st:dim:257-st, st:dim:257-st, real(max_ilp), -imag(max_ilp), 2, 'c'); 
    title(['Level ', num2str(level), ' ICP coefficients']); 
    
end

