%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to test the transform of a dual-tree into ILP and ICP coefficients
% and the feature vectors formed from them

% Using lenna as a test image, compute the ICP and ILP transforms of the
% dual-tree, and display the weighted-sum of the maximal coefficients at
% each level

lenna = u_load('C:\isbe\matlab_code\trunk\wavelet_dtcwt_toolbox4_1\lenna.mat');
dt_lenna = dtwavexfm2(lenna, 5, 'near_sym_b','qshift_b');
%%
[ilp icp] = mb_dual_tree_transform(dt_lenna);
%
weighted_icp = max(icp{1}, [], 3);
weighted_ilp = max(ilp{1}, [], 3);
%imwrite(complex2rgb(weighted_icp, pi), 'C:\isbe\dev\background\dt_transforms\lenna_max_icp_level1.bmp');
imwrite(complex2rgb(weighted_ilp, [-pi/2 pi/2]), 'C:\isbe\dev\background\dt_transforms\lenna_max_ilp_level1.bmp');

%
for level = 2:4
    %Get maximal coefficients across the sub-band across for this level
    max_icp = max(icp{level}, [], 3);
    max_ilp = max(ilp{level}, [], 3);
    
    %imwrite(complex2rgb(max_icp, pi), ['C:\isbe\dev\background\dt_transforms\lenna_max_icp_level', num2str(level), '.bmp']);
    imwrite(complex2rgb(max_ilp, [-pi/2 pi/2]), ['C:\isbe\dev\background\dt_transforms\lenna_max_ilp_level', num2str(level), '.bmp']);
    
    %Generate upscaling indices
    [R C] = size(max_icp);
    idx = kron(reshape(1:R*C,R,C), ones(2^(level-1)));
    
    %Add to the weighted sums
    weighted_icp = weighted_icp + max_icp(idx);
    weighted_ilp = weighted_ilp + max_ilp(idx);
end

%figure; image(complex2rgb(weighted_icp, pi)); axis image;
figure; image(complex2rgb(weighted_ilp, [-pi/2 pi/2])); axis image;
%
%imwrite(complex2rgb(weighted_icp, pi), 'C:\isbe\dev\background\dt_transforms\lenna_max_icp_weighted.bmp');
imwrite(complex2rgb(weighted_ilp, [-pi/2 pi/2]), 'C:\isbe\dev\background\dt_transforms\lenna_max_ilp_weighted.bmp');

%%
% Do the same but for a nice mammographic patch
dt_mammo = u_load('C:\isbe\dev\background\dual_tree\normal512\o04_011RML_dual_tree');

[ilp icp] = mb_dual_tree_transform(dt_mammo);
%
weighted_icp = max(icp{1}, [], 3);
weighted_ilp = max(ilp{1}, [], 3);
%imwrite(complex2rgb(weighted_icp), 'C:\isbe\dev\background\dt_transforms\mammo_max_icp_level1.bmp');
imwrite(complex2rgb(weighted_ilp, [-pi/2 pi/2]), 'C:\isbe\dev\background\dt_transforms\mammo_max_ilp_level1.bmp');

%
for level = 2:4
    %Get maximal coefficients across the sub-band across for this level
    max_icp = max(icp{level}, [], 3);
    max_ilp = max(ilp{level}, [], 3);
    
    %figure; image(complex2rgb(max_icp, pi)); axis image;
    figure; image(complex2rgb(max_ilp, [-pi/2 pi/2])); axis image;
    
    %imwrite(complex2rgb(max_icp, pi), ['C:\isbe\dev\background\dt_transforms\mammo_max_icp_level', num2str(level), '.bmp']);
    imwrite(complex2rgb(max_ilp, [-pi/2 pi/2]), ['C:\isbe\dev\background\dt_transforms\mammo_max_ilp_level', num2str(level), '.bmp']);
    
    %Generate upscaling indices
    [R C] = size(max_icp);
    idx = kron(reshape(1:R*C,R,C), ones(2^(level-1)));
    
    %Add to the weighted sums
    weighted_icp = weighted_icp + max_icp(idx);
    weighted_ilp = weighted_ilp + max_ilp(idx);
end

%figure; image(complex2rgb(weighted_icp, pi)); axis image;
figure; image(complex2rgb(weighted_ilp, [-pi/2 pi/2])); axis image;
%
%imwrite(complex2rgb(weighted_icp, pi), 'C:\isbe\dev\background\dt_transforms\mammo_max_icp_weighted.bmp');
imwrite(complex2rgb(weighted_ilp, [-pi/2 pi/2]), 'C:\isbe\dev\background\dt_transforms\mammo_max_ilp_weighted.bmp');