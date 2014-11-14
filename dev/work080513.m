%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Work May 14th 2008
%
% Visual analysis of Dual-Tree decompostion sub-bands -is there anything
% there? Are Steerable Pyramid subbands really more intuitive to look at -
% or ami deluding myself?
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load image a decompose
I = imread('C:\isbe\dev\background\images\normal512\o04_001LCC_1024_3427_865.bmp');
[p p_sizes] = mb_buildSFpyr(I, 5, 3);
pyr = mb_change_pyramid_form(p, p_sizes);

[Yl,Yh] = dtwavexfm2(double(I),5,'near_sym_b','qshift_b');

%%
% Show magnitude, phase, real and imaginary parts of complex dual-tree
% wavelet coefficients
for level = 1:4
    for ori = 1:6
        figure; 
        subplot(2,2,1); imagesc(abs(Yh{level}(:,:,ori))); axis image; 
        subplot(2,2,2); imagesc(angle(Yh{level}(:,:,ori))); axis image; 
        subplot(2,2,3); imagesc(real(Yh{level}(:,:,ori))); axis image; 
        subplot(2,2,4); imagesc(imag(Yh{level}(:,:,ori))); axis image;
    end
end

%%
% Show coefficients of Steerable Pyramid
for level = 2:6
    for ori = 1:4
        figure; imagesc(pyr{level,ori}); axis image; colormap(jet(256));
    end
end
