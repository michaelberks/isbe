% We've got a method for extracting a feature vector for any point in an
% image. Now test to see if this has the properties we want. i.e. Rotation
% and shift invariance

% One simple method. Take a test image, rotate it through 90, 180 and 270
% degrees (so we know the pixel values remain exact) and compare the
% feature vectors for equivalent points

% Equivalent points for an RxC image are as follows:
% (r, c) -> r90(r, c) = (C-c+1, r)
%        -> r180(r, c) = (R-r+1, C-c+1)
%        -> r270(r, c) = (c, R-r+1)

% load lenna and create rotated copies
lenna = u_load('C:\isbe\matlab_code\trunk\wavelet_dtcwt_toolbox4_1\lenna.mat');
lenna90 = imrotate(lenna, 90);
lenna180 = imrotate(lenna, 180);
lenna270 = imrotate(lenna, 270);
[R C] = size(lenna);

% Show the images if necessary
% figure; imagesc(lenna); axis image; colormap(gray(256));
% figure; imagesc(lenna90); axis image; colormap(gray(256));
% figure; imagesc(lenna180); axis image; colormap(gray(256)); 
% figure; imagesc(lenna270); axis image; colormap(gray(256));

% compute dual-tree for each image
dt_lenna = dtwavexfm2(lenna, 4, 'near_sym_b','qshift_b');
dt_lenna90 = dtwavexfm2(lenna90, 4, 'near_sym_b','qshift_b');
dt_lenna180 = dtwavexfm2(lenna180, 4, 'near_sym_b','qshift_b');
dt_lenna270 = dtwavexfm2(lenna270, 4, 'near_sym_b','qshift_b');
%%
% compute ILP and ICP coefficients for each image
icp_0 = inter_coefficient_product(dt_lenna);
icp_90 = inter_coefficient_product(dt_lenna90);
icp_180 = inter_coefficient_product(dt_lenna180);
icp_270 = inter_coefficient_product(dt_lenna270);
%%
ilp_0 = dtcwt2ilp3(dt_lenna);
ilp_90 = dtcwt2ilp3(dt_lenna90);
ilp_180 = dtcwt2ilp3(dt_lenna180);
ilp_270 = dtcwt2ilp3(dt_lenna270);
%%
for lev = 1:3
    for ori = 1:3
        
        ilp_diff_90{lev}(:,:,ori) = ilp_0{lev}(:,:,ori) - rot90(ilp_90{lev}(:,:,ori+3),-1);
        ilp_diff_180{lev}(:,:,ori) = ilp_0{lev}(:,:,ori) - rot90(ilp_180{lev}(:,:,ori),-2);
        ilp_diff_270{lev}(:,:,ori) = ilp_0{lev}(:,:,ori) - rot90(ilp_270{lev}(:,:,ori+3),-3);
        
    end
    for ori = 4:6
        
        ilp_diff_90{lev}(:,:,ori) = ilp_0{lev}(:,:,ori) - rot90(ilp_90{lev}(:,:,ori-3),-1);
        ilp_diff_180{lev}(:,:,ori) = ilp_0{lev}(:,:,ori) - rot90(ilp_180{lev}(:,:,ori),-2);
        ilp_diff_270{lev}(:,:,ori) = ilp_0{lev}(:,:,ori) - rot90(ilp_270{lev}(:,:,ori-3),-3);
        
    end
end
%%
for lev = 1:3
    f90 = figure;
    f180 = figure;
    f270 = figure;
    for ori = 1:6
        figure(f90); subplot(2,3,ori); imagesc(abs(ilp_diff_90{lev}(:,:,ori))); axis image; colormap(hsv); colorbar;
        figure(f180); subplot(2,3,ori); imagesc(abs(ilp_diff_180{lev}(:,:,ori))); axis image; colormap(hsv); colorbar;
        figure(f270); subplot(2,3,ori); imagesc(abs(ilp_diff_270{lev}(:,:,ori))); axis image; colormap(hsv); colorbar;
    end
end
%%
% Get maximum subband
subband_idx0 = cell(4, 1);
subband_idx90 = cell(4, 1);
subband_idx180 = cell(4, 1);
subband_idx270 = cell(4, 1);
for lev = 1:4
    [dummy subband_idx0{lev}] = max(dt_lenna{lev}, [], 3);
    [dummy subband_idx90{lev}] = max(dt_lenna90{lev}, [], 3);
    [dummy subband_idx180{lev}] = max(dt_lenna180{lev}, [], 3);
    [dummy subband_idx270{lev}] = max(dt_lenna270{lev}, [], 3);
    clear dummy;
end
%%
%check they're symmetric
for lev = 1:4
    sum(sum(mod(subband_idx0{lev} - rot90(subband_idx90{lev}, -1),3)))
    sum(sum(mod(subband_idx0{lev} - rot90(subband_idx180{lev}, -2),6)))
    sum(sum(mod(subband_idx0{lev} - rot90(subband_idx270{lev}, -3),3)))
end
%%
% We pick (82 162) as interesting point - edge of the hat
dtc_l0 = mb_get_dt_vector('pts', [82 162], 'dual_tree', dt_lenna);
dtc_l90 = mb_get_dt_vector('pts', [95 82], 'dual_tree', dt_lenna90);
dtc_l180 = mb_get_dt_vector('pts', [175 95], 'dual_tree', dt_lenna180);
dtc_l270 = mb_get_dt_vector('pts', [162 175], 'dual_tree', dt_lenna270);
%%
% Plot the vectors in a single figure;
figure;
subplot(2,2,1); mb_plot_dt_vector(dtc_l0);
subplot(2,2,2); mb_plot_dt_vector(dtc_l90);
subplot(2,2,3); mb_plot_dt_vector(dtc_l180);
subplot(2,2,4); mb_plot_dt_vector(dtc_l270);
%%
%Plot the feature vectors as colorplots - this method works nicely and is
%suitable for comparing a large numer of points
figure;
subplot(4,1,1); image(complex2rgb( dtc_l0(:,1:3:end-2).*exp(i*dtc_l0(:,2:3:end-1)) ));
subplot(4,1,2); image(complex2rgb( dtc_l90(:,1:3:end-2).*exp(i*dtc_l90(:,2:3:end-1)) ));
subplot(4,1,3); image(complex2rgb( dtc_l180(:,1:3:end-2).*exp(i*dtc_l180(:,2:3:end-1)) ));
subplot(4,1,4); image(complex2rgb( dtc_l270(:,1:3:end-2).*exp(i*dtc_l270(:,2:3:end-1)) ));

%%
% So lets try it on a number of randomly selected points
%[r0 c0] = ind2sub([R C], randsample(1:R*C, 64)');
[r0 c0] = ind2sub([R C], (1:R*C)');

r90 = C-c0+1; c90 = r0;
r180 = R-r0+1; c180 = C-c0+1;
r270 = c0; c270 = R-r0+1;
%%
%profile on;
dtc_l0 = mb_get_dt_vector('pts', [r0 c0], 'dual_tree', dt_lenna);
%profile viewer;
dtc_l90 = mb_get_dt_vector('pts', [r90 c90], 'dual_tree', dt_lenna90);
dtc_l180 = mb_get_dt_vector('pts', [r180 c180], 'dual_tree', dt_lenna180);
dtc_l270 = mb_get_dt_vector('pts', [r270 c270], 'dual_tree', dt_lenna270);
%%
figure;
subplot(1,4,1); image(complex2rgb( dtc_l0(:,1:3:end-2).*exp(i*dtc_l0(:,2:3:end-1)) ));
subplot(1,4,2); image(complex2rgb( dtc_l90(:,1:3:end-2).*exp(i*dtc_l90(:,2:3:end-1)) ));
subplot(1,4,3); image(complex2rgb( dtc_l180(:,1:3:end-2).*exp(i*dtc_l180(:,2:3:end-1)) ));
subplot(1,4,4); image(complex2rgb( dtc_l270(:,1:3:end-2).*exp(i*dtc_l270(:,2:3:end-1)) ));
%%
% Check differences in magnitude
diff_m90 = dtc_l90(:,1:3:end-2) - dtc_l0(:,1:3:end-2);
diff_m180 = dtc_l180(:,1:3:end-2) - dtc_l0(:,1:3:end-2);
diff_m270 = dtc_l270(:,1:3:end-2) - dtc_l0(:,1:3:end-2);

% Check differences in orientation (i.e. from ICP)
diff_o90 = dtc_l90(:,2:3:end-1) - dtc_l0(:,2:3:end-1);
diff_o180 = dtc_l180(:,2:3:end-1) - dtc_l0(:,2:3:end-1);
diff_o270 = dtc_l270(:,2:3:end-1) - dtc_l0(:,2:3:end-1);

%
% Check differences in scale phase (i.e. from ILP)
diff_s90 = dtc_l90(:,3:3:end-3) - dtc_l0(:,3:3:end-3);
diff_s180 = dtc_l180(:,3:3:end-3) - dtc_l0(:,3:3:end-3);
diff_s270 = dtc_l270(:,3:3:end-3) - dtc_l0(:,3:3:end-3);

%%
max_diff = max(max(abs([diff90; diff180; diff270])));
figure;
subplot(1,3,1); image(complex2rgb( diff90, max_diff ));
subplot(1,3,2); image(complex2rgb( diff180, max_diff ));
subplot(1,3,3); image(complex2rgb( diff270, max_diff ));

%%
figure;
subplot(2,2,1); mb_plot_dt_vector(dtc_l0(46,:));
subplot(2,2,2); mb_plot_dt_vector(dtc_l90(46,:));
subplot(2,2,3); mb_plot_dt_vector(dtc_l180(46,:));
subplot(2,2,4); mb_plot_dt_vector(dtc_l270(46,:));
%%
for level = 1:3
    figure;
    subplot(1,3,1); imagesc(reshape(mod(diff_s90(:,level), pi/2), 256, 256)); axis image; colormap(hsv); colorbar;
    subplot(1,3,2); imagesc(reshape(mod(diff_s180(:,level), pi/2), 256, 256)); axis image; colormap(hsv); colorbar
    subplot(1,3,3); imagesc(reshape(mod(diff_s270(:,level), pi/2), 256, 256)); axis image; colormap(hsv); colorbar;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
r0 = 83; c0 = 192;
r90 = C-c0+1; c90 = r0;
r180 = R-r0+1; c180 = C-c0+1;
r270 = c0; c270 = R-r0+1;

for lev = 1:4
    [subband_idx0{lev}(ceil(r0 / 2^lev), ceil(c0 / 2^lev)) ...
     subband_idx90{lev}(ceil(r90 / 2^lev), ceil(c90 / 2^lev)) ...
     subband_idx180{lev}(ceil(r180 / 2^lev), ceil(c180 / 2^lev)) ...
     subband_idx270{lev}(ceil(r270 / 2^lev), ceil(c270 / 2^lev)) ]
end
%%
for lev = 1:3
    [ilp_0{lev}(ceil(r0 / 2^lev), ceil(c0 / 2^lev), subband_idx0{lev}(ceil(r0 / 2^lev), ceil(c0 / 2^lev))) ...
     ilp_90{lev}(ceil(r90 / 2^lev), ceil(c90 / 2^lev), subband_idx90{lev}(ceil(r90 / 2^lev), ceil(c90 / 2^lev))) ...
     ilp_180{lev}(ceil(r180 / 2^lev), ceil(c180 / 2^lev), subband_idx180{lev}(ceil(r180 / 2^lev), ceil(c180 / 2^lev))) ...
     ilp_270{lev}(ceil(r270 / 2^lev), ceil(c270 / 2^lev), subband_idx270{lev}(ceil(r270 / 2^lev), ceil(c270 / 2^lev))) ]
end
%%
for lev = 1:3
    for ori = 1:3
        display(['Level ' num2str(lev), ', subband ' num2str(ori), ': ']);
        [ilp_0{lev}(ceil(r0 / 2^lev), ceil(c0 / 2^lev), ori) ...
         ilp_90{lev}(ceil(r90 / 2^lev), ceil(c90 / 2^lev), ori+3) ...
         ilp_180{lev}(ceil(r180 / 2^lev), ceil(c180 / 2^lev), ori) ...
         ilp_270{lev}(ceil(r270 / 2^lev), ceil(c270 / 2^lev), ori+3) ]
    end
    for ori = 4:6
        display(['Level ' num2str(lev), ', subband ' num2str(ori), ': '])
        [ilp_0{lev}(ceil(r0 / 2^lev), ceil(c0 / 2^lev), ori) ...
         ilp_90{lev}(ceil(r90 / 2^lev), ceil(c90 / 2^lev), ori-3) ...
         ilp_180{lev}(ceil(r180 / 2^lev), ceil(c180 / 2^lev), ori) ...
         ilp_270{lev}(ceil(r270 / 2^lev), ceil(c270 / 2^lev), ori-3) ]
    end
end
%%
for lev = 1:3
    for ori = 1:3
        display(['Level ' num2str(lev), ', subband ' num2str(ori), ': ']);
        [ilp_0{lev}(ceil(r0 / 2^lev), ceil(c0 / 2^lev), ori) ...
         ilp_90{lev}(ceil(r90 / 2^lev), ceil(c90 / 2^lev), ori+3) ...
         ilp_180{lev}(ceil(r180 / 2^lev), ceil(c180 / 2^lev), ori) ...
         ilp_270{lev}(ceil(r270 / 2^lev), ceil(c270 / 2^lev), ori+3) ]
    end
    for ori = 4:6
        display(['Level ' num2str(lev), ', subband ' num2str(ori), ': '])
        [ilp_0{lev}(ceil(r0 / 2^lev), ceil(c0 / 2^lev), ori) ...
         ilp_90{lev}(ceil(r90 / 2^lev), ceil(c90 / 2^lev), ori-3) ...
         ilp_180{lev}(ceil(r180 / 2^lev), ceil(c180 / 2^lev), ori) ...
         ilp_270{lev}(ceil(r270 / 2^lev), ceil(c270 / 2^lev), ori-3) ]
    end
end
%%
for lev = 1:3
    for ori = 1:3
        display(['Level ' num2str(lev), ', subband ' num2str(ori), ': ']);
        [ilp_0{lev}(ceil(r0 / 2^lev), ceil(c0 / 2^lev), ori) ...
         new_ilp90{lev}(ceil(r90 / 2^lev), ceil(c90 / 2^lev), ori+3) ...
         new_ilp180{lev}(ceil(r180 / 2^lev), ceil(c180 / 2^lev), ori) ...
         new_ilp270{lev}(ceil(r270 / 2^lev), ceil(c270 / 2^lev), ori+3) ]
    end
    for ori = 4:6
        display(['Level ' num2str(lev), ', subband ' num2str(ori), ': '])
        [ilp_0{lev}(ceil(r0 / 2^lev), ceil(c0 / 2^lev), ori) ...
         new_ilp90{lev}(ceil(r90 / 2^lev), ceil(c90 / 2^lev), ori-3) ...
         new_ilp180{lev}(ceil(r180 / 2^lev), ceil(c180 / 2^lev), ori) ...
         new_ilp270{lev}(ceil(r270 / 2^lev), ceil(c270 / 2^lev), ori-3) ]
    end
end
%%
for lev = 1:3
    for ori = 1:3
        new_ilp90{lev}(:,:,ori) = complex(-imag(ilp_90{lev}(:,:,ori)), -real(ilp_90{lev}(:,:,ori)));
        new_ilp180{lev}(:,:,ori) = complex(-real(ilp_180{lev}(:,:,ori)), imag(ilp_180{lev}(:,:,ori)));
        new_ilp270{lev}(:,:,ori) = complex(-imag(ilp_270{lev}(:,:,ori)), real(ilp_270{lev}(:,:,ori)));
    end
    for ori = 4:6
        new_ilp90{lev}(:,:,ori) = complex(imag(ilp_90{lev}(:,:,ori)), -real(ilp_90{lev}(:,:,ori)));
        new_ilp180{lev}(:,:,ori) = complex(real(ilp_180{lev}(:,:,ori)), -imag(ilp_180{lev}(:,:,ori)));
        new_ilp270{lev}(:,:,ori) = complex(-imag(ilp_270{lev}(:,:,ori)), -real(ilp_270{lev}(:,:,ori)));
    end
end
%%
for lev = 1:3
    for ori = 1:3
        if(sum(sum(abs(ilp_0{lev}(:,:,ori) - rot90(new_ilp90{lev}(:,:,ori+3), -1)))) > 1e-4); display(['90: Level ' num2str(lev), ', subband ' num2str(ori)]); end
        if(sum(sum(abs(ilp_0{lev}(:,:,ori) - rot90(new_ilp180{lev}(:,:,ori), -2)))) > 1e-4); display(['180: Level ' num2str(lev), ', subband ' num2str(ori)]); end
        if(sum(sum(abs(ilp_0{lev}(:,:,ori) - rot90(new_ilp270{lev}(:,:,ori+3), -3)))) > 1e-4); display(['270: Level ' num2str(lev), ', subband ' num2str(ori)]); end
    end
    for ori = 4:6
        if(sum(sum(abs(ilp_0{lev}(:,:,ori) - rot90(new_ilp90{lev}(:,:,ori-3), -1)))) > 1e-4); display(['90: Level ' num2str(lev), ', subband ' num2str(ori)]); end
        if(sum(sum(abs(ilp_0{lev}(:,:,ori) - rot90(new_ilp180{lev}(:,:,ori), -2)))) > 1e-4); display(['180: Level ' num2str(lev), ', subband ' num2str(ori)]); end
        if(sum(sum(abs(ilp_0{lev}(:,:,ori) - rot90(new_ilp270{lev}(:,:,ori-3), -3)))) > 1e-4); display(['270: Level ' num2str(lev), ', subband ' num2str(ori)]); end
    end
end
%%
w = [-3 -1; -3 -3; -1 -3; 1 -3; 3 -3; 3 -1]*pi/2.15; % Nominally j * pi/2, but reduced a bit due to asymmetry of subband freq responses.
w2 = [-3 -1; -3 -3; -1 -3; -1 3; -3 3; -3 1]*pi/2.15;
p = [1 3]/4;

for lev = 1:3
    for ori = 1:6
        dt_int0{lev}(:,:,ori) = cpxinterp2(dt_lenna{lev+1}(:,:,ori), p-0.5, w(ori,:),'spline');
        dt_int90{lev}(:,:,ori) = cpxinterp2(dt_lenna90{lev+1}(:,:,ori), p-0.5, w(ori,:),'spline');
        dt_int180{lev}(:,:,ori) = cpxinterp2(dt_lenna180{lev+1}(:,:,ori), p-0.5, w(ori,:),'spline');
        dt_int270{lev}(:,:,ori) = cpxinterp2(dt_lenna270{lev+1}(:,:,ori), p-0.5, w(ori,:),'spline');
        
        new_dt_int0{lev}(:,:,ori) = cpxinterp2(new_dt_0{lev+1}(:,:,ori), p-0.5, w(ori,:),'spline');
        new_dt_int90{lev}(:,:,ori) = cpxinterp2(new_dt_90{lev+1}(:,:,ori), p-0.5, w(ori,:),'spline');
        new_dt_int180{lev}(:,:,ori) = cpxinterp2(new_dt_180{lev+1}(:,:,ori), p-0.5, w(ori,:),'spline');
        new_dt_int270{lev}(:,:,ori) = cpxinterp2(new_dt_270{lev+1}(:,:,ori), p-0.5, w(ori,:),'spline');
        
        new_dt_int0b{lev}(:,:,ori) = cpxinterp2(new_dt_0{lev+1}(:,:,ori), p-0.5, w2(ori,:),'spline');
        new_dt_int90b{lev}(:,:,ori) = cpxinterp2(new_dt_90{lev+1}(:,:,ori), p-0.5, w2(ori,:),'spline');
        new_dt_int180b{lev}(:,:,ori) = cpxinterp2(new_dt_180{lev+1}(:,:,ori), p-0.5, w2(ori,:),'spline');
        new_dt_int270b{lev}(:,:,ori) = cpxinterp2(new_dt_270{lev+1}(:,:,ori), p-0.5, w2(ori,:),'spline');
    end
    
    for ori = 1:3
       
        new_dt_int0a{lev}(:,:,ori) = dt_int0{lev}(:,:,ori);
        new_dt_int90a{lev}(:,:,ori) = dt_int90{lev}(:,:,ori);
        new_dt_int180a{lev}(:,:,ori) = dt_int180{lev}(:,:,ori);
        new_dt_int270a{lev}(:,:,ori) = dt_int270{lev}(:,:,ori);
    end
    for ori = 4:6
        new_dt_int0a{lev}(:,:,ori) = complex(imag(dt_int0{lev}(:,:,ori)), real(dt_int0{lev}(:,:,ori)));
        new_dt_int90a{lev}(:,:,ori) = complex(imag(dt_int90{lev}(:,:,ori)), real(dt_int90{lev}(:,:,ori)));
        new_dt_int180a{lev}(:,:,ori) = complex(imag(dt_int180{lev}(:,:,ori)), real(dt_int180{lev}(:,:,ori)));
        new_dt_int270a{lev}(:,:,ori) = complex(imag(dt_int270{lev}(:,:,ori)), real(dt_int270{lev}(:,:,ori)));
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
new_dt_0 = dt_lenna;
new_dt_90 = dt_lenna90;
new_dt_180 = dt_lenna180;
new_dt_270 = dt_lenna270;

for lev = 1:4
    for ori = 4:6
        new_dt_0{lev}(:,:,ori) = complex(imag(dt_lenna{lev}(:,:,ori)), real(dt_lenna{lev}(:,:,ori)));
        new_dt_90{lev}(:,:,ori) = complex(imag(dt_lenna90{lev}(:,:,ori)), real(dt_lenna90{lev}(:,:,ori)));
        new_dt_180{lev}(:,:,ori) = complex(imag(dt_lenna180{lev}(:,:,ori)), real(dt_lenna180{lev}(:,:,ori)));
        new_dt_270{lev}(:,:,ori) = complex(imag(dt_lenna270{lev}(:,:,ori)), real(dt_lenna270{lev}(:,:,ori)));
    end
end
%%
[reshape(dt_lenna{3}(1,1,:), 1,[]); reshape(dt_lenna90{3}(32,1,[4:6, 1:3]), 1,[]); reshape(dt_lenna180{3}(32,32,:), 1,[]); reshape(dt_lenna270{3}(1,32,[4:6, 1:3]), 1,[])]
[reshape(new_dt_0{3}(1,1,:), 1,[]); reshape(new_dt_90{3}(32,1,[4:6, 1:3]), 1,[]); reshape(new_dt_180{3}(32,32,:), 1,[]); reshape(new_dt_270{3}(1,32,[4:6, 1:3]), 1,[])]
%%
angle([reshape(dt_lenna{3}(1,1,:), 1,[]); reshape(dt_lenna90{3}(32,1,[4:6, 1:3]), 1,[]); reshape(dt_lenna180{3}(32,32,:), 1,[]); reshape(dt_lenna270{3}(1,32,[4:6, 1:3]), 1,[])])
angle([reshape(new_dt_0{3}(1,1,:), 1,[]); reshape(new_dt_90{3}(32,1,[4:6, 1:3]), 1,[]); reshape(new_dt_180{3}(32,32,:), 1,[]); reshape(new_dt_270{3}(1,32,[4:6, 1:3]), 1,[])])

angle([reshape(dt_lenna{3}(1,1,:), 1,[]); reshape(dt_lenna90{3}(32,1,[4:6, 1:3]), 1,[]); reshape(dt_lenna180{3}(32,32,:), 1,[]); reshape(dt_lenna270{3}(1,32,[4:6, 1:3]), 1,[])])
angle([reshape(new_dt_0{3}(1,1,:), 1,[]); reshape(new_dt_90{3}(32,1,[4:6, 1:3]), 1,[]); reshape(new_dt_180{3}(32,32,:), 1,[]); reshape(new_dt_270{3}(1,32,[4:6, 1:3]), 1,[])])

angle([reshape(dt_lenna{3}(1,1,:), 1,[]); reshape(dt_lenna90{3}(32,1,[4:6, 1:3]), 1,[]); reshape(dt_lenna180{3}(32,32,:), 1,[]); reshape(dt_lenna270{3}(1,32,[4:6, 1:3]), 1,[])])
angle([reshape(new_dt_0{3}(1,1,:), 1,[]); reshape(new_dt_90{3}(32,1,[4:6, 1:3]), 1,[]); reshape(new_dt_180{3}(32,32,:), 1,[]); reshape(new_dt_270{3}(1,32,[4:6, 1:3]), 1,[])])

angle([reshape(dt_lenna{3}(1,1,:), 1,[]); reshape(dt_lenna90{3}(32,1,[4:6, 1:3]), 1,[]); reshape(dt_lenna180{3}(32,32,:), 1,[]); reshape(dt_lenna270{3}(1,32,[4:6, 1:3]), 1,[])])
angle([reshape(new_dt_0{3}(1,1,:), 1,[]); reshape(new_dt_90{3}(32,1,[4:6, 1:3]), 1,[]); reshape(new_dt_180{3}(32,32,:), 1,[]); reshape(new_dt_270{3}(1,32,[4:6, 1:3]), 1,[])])

angle([reshape(dt_lenna{3}(1,1,:), 1,[]); reshape(dt_lenna90{3}(32,1,[4:6, 1:3]), 1,[]); reshape(dt_lenna180{3}(32,32,:), 1,[]); reshape(dt_lenna270{3}(1,32,[4:6, 1:3]), 1,[])])
angle([reshape(new_dt_0{3}(1,1,:), 1,[]); reshape(new_dt_90{3}(32,1,[4:6, 1:3]), 1,[]); reshape(new_dt_180{3}(32,32,:), 1,[]); reshape(new_dt_270{3}(1,32,[4:6, 1:3]), 1,[])])

angle([reshape(dt_lenna{3}(1,1,:), 1,[]); reshape(dt_lenna90{3}(32,1,[4:6, 1:3]), 1,[]); reshape(dt_lenna180{3}(32,32,:), 1,[]); reshape(dt_lenna270{3}(1,32,[4:6, 1:3]), 1,[])])
angle([reshape(new_dt_0{3}(1,1,:), 1,[]); reshape(new_dt_90{3}(32,1,[4:6, 1:3]), 1,[]); reshape(new_dt_180{3}(32,32,:), 1,[]); reshape(new_dt_270{3}(1,32,[4:6, 1:3]), 1,[])])
%%
new_icp_0 = inter_coefficient_product4(new_dt_0);
new_icp_90 = inter_coefficient_product4(new_dt_90);
new_icp_180 = inter_coefficient_product4(new_dt_180);
new_icp_270 = inter_coefficient_product4(new_dt_270);
%%
for level = 1:4
    dim = 2^level;
    st = (dim+1) / 2;
    [max_icp] = max(new_icp_0{level}, [], 3);
    figure; imagesc(lenna); axis image; colormap(gray(256)); hold on;
    quiver(st:dim:256-st+1, st:dim:256-st+1, real(max_icp), -imag(max_icp));
end
%%
for lev = 1:4
    for ori = 1:3
        temp_diff = mod(rot90(angle(new_icp_90{lev}(:,:,ori+3)),-1) - pi/2 - angle(new_icp_0{lev}(:,:,ori)), pi);
        temp_diff(temp_diff < 1e-4) = 0;
        figure; imagesc(temp_diff); axis image; colormap(hsv); colorbar;
    end
    for ori = 4:6
        temp_diff = mod(rot90(angle(new_icp_90{lev}(:,:,ori-3)),-1) - pi/2 - angle(new_icp_0{lev}(:,:,ori)), pi);
        temp_diff(temp_diff < 1e-4) = 0;
        figure; imagesc(temp_diff); axis image; colormap(hsv); colorbar;
    end
end
%%
for lev = 1:4
    for ori = 1:6
        temp_diff = mod(rot90(angle(new_icp_180{lev}(:,:,ori)),-2) - angle(new_icp_0{lev}(:,:,ori)), pi);
        temp_diff(temp_diff < 1e-4) = 0;
        figure; imagesc(temp_diff); axis image; colormap(hsv); colorbar;
    end
end
%%
for lev = 1:4
    for ori = 1:3
        temp_diff = mod(rot90(angle(new_icp_270{lev}(:,:,ori+3)),-3) + pi/2 - angle(new_icp_0{lev}(:,:,ori)), pi);
        temp_diff(temp_diff < 1e-4) = 0;
        figure; imagesc(temp_diff); axis image; colormap(hsv);
    end
    for ori = 4:6
        temp_diff = mod(rot90(angle(new_icp_270{lev}(:,:,ori-3)),-3) + pi/2 - angle(new_icp_0{lev}(:,:,ori)), pi);
        temp_diff(temp_diff < 1e-4) = 0;
        figure; imagesc(temp_diff); axis image; colormap(hsv); colorbar;
    end
end
%%
for lev = 1:4
    max_icp0 = max(new_icp_0{lev}, [], 3);
    max_icp90 = max(new_icp_90{lev}, [], 3);
    max_icp180 = max(new_icp_180{lev}, [], 3);
    max_icp270 = max(new_icp_270{lev}, [], 3);
    
    new_diff_o90 = pi_rotate(rot90(angle(max_icp90),-1) - pi/2 - angle(max_icp0));
    new_diff_o180 = pi_rotate(rot90(angle(max_icp180),-2) - angle(max_icp0));
    new_diff_o270 = pi_rotate(rot90(angle(max_icp270),-3) + pi/2 - angle(max_icp0));
    
    figure; imagesc(new_diff_o90); axis image; colormap(hsv); caxis([0 2*pi]);
    figure; imagesc(new_diff_o180); axis image; colormap(hsv); caxis([0 2*pi]);
    figure; imagesc(new_diff_o270); axis image; colormap(hsv); caxis([0 2*pi]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% So we're happy the orientation calculation is now sorted - what about
% inter-scale phase
%%
new_ilp_0 = dtcwt2ilp(new_dt_0);
new_ilp_90 = dtcwt2ilp(new_dt_90);
new_ilp_180 = dtcwt2ilp(new_dt_180);
new_ilp_270 = dtcwt2ilp(new_dt_270);
%%
angle([reshape(new_ilp_0{3}(1,1,:), 1,[]);...
    reshape(new_ilp_90{3}(32,1,[4:6, 1:3]), 1,[]);...
    reshape(new_ilp_180{3}(32,32,:), 1,[]);...
    reshape(new_ilp_270{3}(1,32,[4:6, 1:3]), 1,[])])
%%
angle([reshape(dt_int0{3}(1,1,:), 1,[]);...
    reshape(dt_int90{3}(32,1,[4:6, 1:3]), 1,[]);...
    reshape(dt_int180{3}(32,32,:), 1,[]);...
    reshape(dt_int270{3}(1,32,[4:6, 1:3]), 1,[])])

angle([reshape(new_dt_int0{3}(1,1,:), 1,[]);...
    reshape(new_dt_int90{3}(32,1,[4:6, 1:3]), 1,[]);...
    reshape(new_dt_int180{3}(32,32,:), 1,[]);...
    reshape(new_dt_int270{3}(1,32,[4:6, 1:3]), 1,[])])
%%
angle([reshape(new_dt_int0a{3}(1,1,:), 1,[]);...
    reshape(new_dt_int90a{3}(32,1,[4:6, 1:3]), 1,[]);...
    reshape(new_dt_int180a{3}(32,32,:), 1,[]);...
    reshape(new_dt_int270a{3}(1,32,[4:6, 1:3]), 1,[])])

angle([reshape(new_dt_int0b{3}(1,1,:), 1,[]);...
    reshape(new_dt_int90b{3}(32,1,[4:6, 1:3]), 1,[]);...
    reshape(new_dt_int180b{3}(32,32,:), 1,[]);...
    reshape(new_dt_int270b{3}(1,32,[4:6, 1:3]), 1,[])])