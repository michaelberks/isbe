mammo = imread('C:\isbe\density\mammograms\018LCC1824.tif');
mammo_small = imresize(mammo, 0.5); clear mammo;
mammo_small = rot90(mammo_small); 
mammo_small = double(mammo_small);
%
figure; imagesc(mammo_small); axis image; colormap(gray(256));
%
dt_mammo = dtwavexfm2(mammo_small, 6);
[ilp_mammo icp_mammo] = mb_dual_tree_transform(dt_mammo, 2);
%%
for level = 2:5
    dt_max{level} = max(dt_mammo{level}, [], 3);
    ilp_max{level} = max(ilp_mammo{level}, [], 3);
    icp_max{level} = max(icp_mammo{level}, [], 3);
    
    dt_mag = histeq(abs(dt_max{level}) / max(abs(dt_max{level}(:))), 1 ./ (1:256));
    ilp_mag = histeq(abs(icp_max{level}) / max(abs(icp_max{level}(:))), 1 ./ (1:256));
    ilp_mag(abs(angle(ilp_max{level})) > 0.4) = 0; 
    icp_mag = histeq(abs(icp_max{level}) / max(abs(icp_max{level}(:))), 1 ./ (1:256));
    figure; 
    imagesc(dt_mag); axis image; colormap(gray(256));
    figure;
    image(complex2rgb(icp_mag.*exp(i*2*angle(icp_max{level})))); axis image;
    figure;
    image(complex2rgb(ilp_mag.*exp(i*angle(ilp_max{level})))); axis image;
end
%%
%Display a figure that shows how far away from zero-phase (and thus an edge) each
%pixel is
figure; imagesc(-abs(angle(ilp_max{2}))); axis image; colormap(gray(256));