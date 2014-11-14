%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Work script for paper: "Determining Multiscale Image Feature Angles from
% Complex Wavelet Phases", Ryan Anderson, Nick Kingsbury, Julien Fauqueur
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Build oriented stepdedges at angles 0:pi/20:19pi/20;
angles = 0:pi/20:19*pi/20;
step_edge_256 = [ones(128,256); 255*ones(128,256)];

step_edges_128 = zeros(128,128,length(angles));

for ii = 1:length(angles);
    
    temp = imrotate(step_edge_256, 180*angles(ii)/pi, 'crop');
    step_edges_128(:,:,ii) = temp(65:192, 65:192);
    clear temp;
    figure; image(step_edges_128(:,:,ii)); colormap(gray(256)); axis image;
end
%%

%First lets play with the pi/10 edge - index 3;

dt_edge_pi10 = dtwavexfm2(step_edges_128(:,:,3), 4, 'near_sym_b','qshift_b');

%%
for level = 1:4
    dim = 2^(6-level);
    axis_lims = [1 dim 1 dim min(abs(dt_edge_pi10{level}(:))) max(abs(dt_edge_pi10{level}(:)))];
    figure;
    for ori = 1:6

        ho(level, ori) = subplot(2,3,ori); mesh(abs(dt_edge_pi10{level}(:,:,ori))); axis(axis_lims); %#ok

    end
    ho_links(level) = linkprop(ho(level,:), {'CameraPosition','CameraUpVector'});%#ok
    
end
%%
icp_edge = inter_coefficient_product2(dt_edge_pi10);
%%

for edge_angle = 1:20
    dt_edge = dtwavexfm2(step_edges_128(:,:,edge_angle), 4, 'near_sym_b','qshift_b');
    dt_edge3 = dt_edge;

    for lev = 1:4
        for ori = 4:6
            dt_edge3{lev}(:,:,ori) = complex(imag(dt_edge{lev}(:,:,ori)), real(dt_edge{lev}(:,:,ori)));
        end
    end
    
    icp_edge3 = inter_coefficient_product4(dt_edge3);
    icp_edge = inter_coefficient_product(dt_edge);
    figure;
    for level = 1:4
        dim = 2^level;
%         figure;
%         for ori = 1:6
%             
%             subplot(2,3,ori);
%             imagesc(step_edges_128(:,:,angle)); axis image; colormap(gray(256)); hold on;
%             quiver(dim:dim:128, dim:dim:128, real(icp_edge{level}(:,:,ori)), -imag(icp_edge{level}(:,:,ori)));
%         end
        [max_icp] = max(icp_edge{level}, [], 3);
        [max_icp3] = max(icp_edge3{level}, [], 3);
        subplot(2,2,level); imagesc(step_edges_128(:,:,edge_angle)); axis image; colormap(gray(256)); hold on;
        quiver(dim:dim:128, dim:dim:128, real(max_icp), -imag(max_icp));
        quiver(dim:dim:128, dim:dim:128, real(max_icp3), -imag(max_icp3), 'r');
    end
end
%%
% Now apply to a set of oriented test lines
test_lines = ones(256);
for ii = 1:8; test_lines((ii*16)-ii+65:(ii*16)+64,:) = 255; end

test_lines_128 = zeros(128,128,length(angles));

figure; imagesc(test_lines); axis image; colormap(gray(256));
for ii = 1:length(angles);
    
    temp = imrotate(test_lines, 180*angles(ii)/pi, 'crop');
    test_lines_128(:,:,ii) = temp(65:192, 65:192);
    clear temp;
    %figure; image(test_lines_128(:,:,ii)); colormap(gray(256)); axis image;
end
%%
for angles = 1:20
    dt_line = dtwavexfm2(test_lines_128(:,:,angles), 4, 'near_sym_b','qshift_b');
    icp_line = inter_coefficient_product(dt_line);
    
    figure;
    for level = 1:4
        dim = 2^level;
        [max_icp indices] = max(icp_line{level}, [], 3);
        
        %subplot(2,2,level); imagesc(test_lines_128(:,:,angle)); axis image; colormap(gray(256)); hold on;
        %quiver(dim:dim:128, dim:dim:128, real(max_icp), -imag(max_icp));
        
        
        nms_icp = non_maximal_supp(abs(max_icp), angle(max_icp));
        subplot(2,2,level); 
        %imagesc(nms_icp>0); axis image; colormap(gray(256)); hold on;
        imagesc(abs(max_icp)); axis image; colormap(gray(256)); hold on;
        max_icp(~nms_icp) = 0;
        quiver(real(max_icp), -imag(max_icp), 1.5);
    end
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% try applying to Lenna test image
load C:\isbe\matlab_code\trunk\wavelet_dtcwt_toolbox4_1\lenna;
dX = dtwavexfm2(X, 4, 'near_sym_b','qshift_b');

icpX = inter_coefficient_product(dX);
%
for level = 1:4
    dim = 2^level;
    
    [max_icp] = max(icpX{level}, [], 3);
    figure; imagesc(X); axis image; colormap(gray(256)); hold on;
    quiver(dim:dim:256, dim:dim:256, real(max_icp), -imag(max_icp), 1.5);
end

%%
%%%%%%%%%%%%%%%%%%%%%%
% Try applying to a mammogram patch
bg_files = dir('C:\isbe\dev\background\images\normal512\*.bmp');
bg = double(imread(['C:\isbe\dev\background\images\normal512\', bg_files(2).name]));

d_bg = dtwavexfm2(bg, 6, 'near_sym_b', 'qshift_b');
icp_bg = inter_coefficient_product(d_bg);
for level = 2:5
    dim = 2^level;
    
    [max_icp] = max(icp_bg{level}, [], 3);
    figure; imagesc(bg); axis image; colormap(gray(256)); hold on;
    quiver(dim:dim:512, dim:dim:512, real(max_icp), -imag(max_icp), 1.5);
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Can we use NMS to select features - probably not yet...
for level = 2:5
    dim = 2^level;
    
    [max_icp] = max(icp_bg{level}, [], 3);
    nms_icp = non_maximal_supp(abs(max_icp), angle(max_icp));
    
    figure; imagesc(nms_icp>0); axis image; colormap(gray(256)); hold on;
    max_icp(~nms_icp) = 0;
    quiver(real(max_icp), -imag(max_icp), 1.5);
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Now try testing the Inter Level Product as defined in the paper: " Coarse
% object recognition using InterLevel Products of Complex Wavelets", also
% by Ryan Anderson, Nick Kingsbury and Julien Fauqueur
%
%

%Build the step_edges above

%First lets look at the ILP tree for 1 step_edge, say the pi/10 angle
dt_edge_pi10 = dtwavexfm2(step_edges_128(:,:,3), 4, 'near_sym_b','qshift_b');
ilp_edge_pi10 = inter_level_product(dt_edge_pi10);

for level = 1:3
    dim = 2^level;
    
    figure;
    for ori = 1:6
        subplot(2,3,ori);
        imagesc(step_edges_128(:,:,3)); axis image; colormap(gray(256)); hold on;
        quiver(dim:dim:128, dim:dim:128, real(ilp_edge_pi10{level}(:,:,ori)), -imag(ilp_edge_pi10{level}(:,:,ori)), 'r');
    end
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do the 1D stuff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Fig 1: Build 1D step-edges
step_edges1D = 2*ones(64, 64);
for offset = 1:64;
    step_edges1D(1:offset, offset) = 1;
end
%
%Construct 1D dual-tree for each offset step-edge
[dummy dt] = dtwavexfm(step_edges1D, 4, 'near_sym_b','qshift_b');
%
%Look at phase of level 1 coefficients - should be linear in support range
% and zero at edge
figure;
for coeff = 1:8
    subplot(4, 2, coeff); hold on; xlabel(['Coeff: ', num2str(coeff)]);
    plot(1:64, angle(dt{1}(coeff,:)), 'b');
    plot(1:64, angle(dt{1}(coeff,:)), 'rx');
    plot([1 64], [0 0], 'r-');
end
%%
figure;
for coeff = 1:8
    subplot(4, 2, coeff); hold on; xlabel(['Coeff: ', num2str(coeff)]);
    plot(1:64, abs(dt{1}(coeff,:)), 'b');
    plot(1:64, abs(dt{1}(coeff,:)), 'rx');
    plot([1 64], [max(abs(dt{1}(coeff,:))) max(abs(dt{1}(coeff,:)))], 'r-');
end

%%
%Look at phase of level 2 coefficients - should be linear in support range
% and zero at edge
figure;
for coeff = 1:8
    subplot(4, 2, coeff); hold on; xlabel(['Coeff: ', num2str(coeff)]);
    plot(1:64, angle(dt{2}(coeff,:)), 'b');
    plot(1:64, angle(dt{2}(coeff,:)), 'rx');
    plot([1 64], [0 0], 'r-');
end
%%
figure;
for coeff = 1:8
    subplot(4, 2, coeff); hold on; xlabel(['Coeff: ', num2str(coeff)]);
    plot(1:64, abs(dt{2}(coeff,:)), 'b');
    plot(1:64, abs(dt{2}(coeff,:)), 'rx');
    plot([1 64], [max(abs(dt{2}(coeff,:))) max(abs(dt{2}(coeff,:)))], 'r-');
end

%%
%Look at phase of level 3 coefficients - should be linear in support range
% and zero at edge
figure;
for coeff = 1:8
    subplot(4, 2, coeff); hold on; xlabel(['Coeff: ', num2str(coeff)]);
    plot(1:64, angle(dt{3}(coeff,:)), 'b');
    plot(1:64, angle(dt{3}(coeff,:)), 'rx');
    plot([1 64], [0 0], 'r-');
end
%%
%Look at mag of level 3 coefficients - should be maximum at edge
figure;
for coeff = 1:8
    subplot(4, 2, coeff); hold on; xlabel(['Coeff: ', num2str(coeff)]);
    plot(1:64, abs(dt{3}(coeff,:)), 'b');
    plot(1:64, abs(dt{3}(coeff,:)), 'rx');
    plot([1 64], [max(abs(dt{3}(coeff,:))) max(abs(dt{3}(coeff,:)))], 'r-');
end
%%
%Look at phase of level 4 coefficients - should be linear in support range
% and zero at edge
figure;
for coeff = 1:4
    subplot(4, 1, coeff); hold on; xlabel(['Coeff: ', num2str(coeff)]);
    plot(1:64, angle(dt{4}(coeff,:)), 'b');
    plot(1:64, angle(dt{4}(coeff,:)), 'rx');
    plot([1 64], [0 0], 'r-');
end
%%
%Look at mag of level 4 coefficients - should be maximum at edge
figure;
for coeff = 1:4
    subplot(4, 1, coeff); hold on; xlabel(['Coeff: ', num2str(coeff)]);
    plot(1:64, abs(dt{4}(coeff,:)), 'b');
    plot(1:64, abs(dt{4}(coeff,:)), 'rx');
    plot([1 64], [max(abs(dt{4}(coeff,:))) max(abs(dt{4}(coeff,:)))], 'r-');
end
%%
[dummy ind1] = max(abs(dt{1}), [], 2);
[dummy ind2] = max(abs(dt{2}), [], 2);
[dummy ind3] = max(abs(dt{3}), [], 2);
[dummy ind4] = max(abs(dt{4}), [], 2);
figure; hold on;
plot(ind1, 4, 'rx', 'markersize', 6.0)
plot(ind2, 3, 'bx', 'markersize', 6.0)
plot(ind3, 2, 'rx', 'markersize', 6.0)
plot(ind4, 1, 'bx', 'markersize', 6.0)
axis([0 65 -2 7]);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fig2: Build the 1D step_edges, I think the section used in the paper would be
%equivalent to
positive_1d_edges = 2*ones(64, 16);
for offset = 1:16;
    positive_1d_edges(1:23+offset, offset) = 1;
end
 
%Construct 1D dual-tree for each offset step-edge
[dummy dt] = dtwavexfm(positive_1d_edges, 4, 'near_sym_b','qshift_b');
%%
negative_1d_edges = ones(64, 16);
for offset = 1:16;
    negative_1d_edges(1:23+offset, offset) = 2;
end
 
%Construct 1D dual-tree for each offset step-edge
[dummy dt] = dtwavexfm(negative_1d_edges, 4, 'near_sym_b','qshift_b');
%%
negative_1d_pulse = 2*ones(64, 16);
for offset = 1:16;
    negative_1d_pulse(23+offset, offset) = 1;
end
 
%Construct 1D dual-tree for each offset step-edge
[dummy dt] = dtwavexfm(negative_1d_pulse, 4, 'near_sym_b','qshift_b');

%%
positive_1d_pulse = ones(64, 16);
for offset = 1:16;
    positive_1d_pulse(23+offset, offset) = 2;
end
 
%Construct 1D dual-tree for each offset step-edge
[dummy dt] = dtwavexfm(positive_1d_pulse, 4, 'near_sym_b','qshift_b');

%%
% Phase double the 4th level co-efficients
dt4_2 = abs(dt{4}).*complex(cos(2*angle(dt{4})), sin(2*angle(dt{4})));

%dim = size(dt{4}, 1);

%Interpolate the 4th level co-efficients
dt4_i = zeros(size(dt{3}, 1), size(dt{3}, 2)); %#ok

% dt_mag = interp1(1:dim, abs(dt{4}), 0.75:0.5:dim+0.25, 'linear', 'extrap');
% dt_phase = interp1(1:dim, angle(dt{4}), 0.75:0.5:dim+0.25, 'linear', 'extrap');
% dt4_i = dt_mag .* exp(i*dt_phase); %#ok
dt4_i = cpxinterp(dt{4}, [-.25 .25], -4.49);

%Phase double the interpolated 4th level co-efficients
dt4_i2 =  abs(dt4_i).*complex(cos(2*angle(dt4_i)), sin(2*angle(dt4_i)));

%Interpolate the phase-doubled 4th level co-efficients

% dt_mag = interp1(1:dim, abs(dt4_2), 0.75:0.5:dim+0.25, 'linear', 'extrap');
% dt_phase = interp1(1:dim, angle(dt4_2), 0.75:0.5:dim+0.25, 'linear', 'extrap');    
% dt4_2i = dt_mag .* exp(i*dt_phase);  %#ok
dt4_2i = cpxinterp(dt4_2, [-.25 .25], -4.4817);

%Compute the ILP coefficients
ilp34_2i = dt{3}.*conj(dt4_2i);
ilp34_i2 = dt{3}.*conj(dt4_i2);
%%
%Display DT coefficients as colorplots
figure; 
subplot(2,2,1); subimage(255*(abs(dt{3})/max(abs(dt{3}(:)))), colormap(gray(256))); axis xy;
title('Level 3 co-efficients: Magnitude');
subplot(2,2,2); subimage(255*(mod(angle(dt{3}),2*pi)/(2*pi)), colormap(hsv(256))); axis xy;
title('Level 3 co-efficients: Phase');
subplot(2,2,3); subimage(complex2rgb(dt{3})); axis xy;
title('Level 3 co-efficients: Magnitude and phase');
subplot(2,2,4); quiver(real(dt{3}), imag(dt{3})); set(gca, 'ytick', 1:8, 'yticklabel', -16:8:40);
title('Level 3 co-efficients: Magnitude and phase');
%
figure; 
subplot(2,2,1); subimage(255*(abs(dt{4})/max(abs(dt{4}(:)))), colormap(gray(256))); axis xy;
title('Level 4 co-efficients: Magnitude');
subplot(2,2,2); subimage(255*(mod(angle(dt{4}),2*pi)/(2*pi)), colormap(hsv(256))); axis xy;
title('Level 4 co-efficients: Phase');
subplot(2,2,3); subimage(complex2rgb(dt{4})); axis xy;
title('Level 4 co-efficients: Magnitude and phase');
subplot(2,2,4); quiver(real(dt{4}), imag(dt{4})); set(gca, 'ytick', 1:8, 'yticklabel', -16:8:40);
title('Level 4 co-efficients: Magnitude and phase');
%
figure; 
subplot(2,2,1); subimage(255*(abs(dt4_i)/max(abs(dt4_i(:)))), colormap(gray(256))); axis xy;
title('Interpolated level 4 co-efficients: Magnitude');
subplot(2,2,2); subimage(255*(mod(angle(dt{3}),2*pi)/(2*pi)), colormap(hsv(256))); axis xy;
title('Interpolated level 4 co-efficients: Phase');
subplot(2,2,3); subimage(complex2rgb(dt4_i)); axis xy;
title('Interpolated level 4 co-efficients: Magnitude and phase');
subplot(2,2,4); quiver(real(dt4_i), imag(dt4_i)); set(gca, 'ytick', 1:8, 'yticklabel', -16:8:40);
title('Interpolated level 4 co-efficients: Magnitude and phase');

figure; 
subplot(2,2,1); subimage(255*(abs(dt4_i2)/max(abs(dt4_i2(:)))), colormap(gray(256))); axis xy;
title('Interpolated then phase-doubled level 4 co-efficients: Magnitude');
subplot(2,2,2); subimage(255*(mod(angle(dt4_i2),2*pi)/(2*pi)), colormap(hsv(256))); axis xy;
title('Interpolated then phase-doubled level 4 co-efficients: Phase');
subplot(2,2,3); subimage(complex2rgb(dt4_i2)); axis xy;
title('Interpolated  then phase-doubled level 4 co-efficients: Magnitude and phase');
subplot(2,2,4); quiver(real(dt4_i2), imag(dt4_i2)); set(gca, 'ytick', 1:8, 'yticklabel', -16:8:40);
title('Interpolated then phase-doubled level 4 co-efficients: Magnitude and phase');

% figure; 
% subplot(2,2,1); subimage(255*(abs(dt4_2i)/max(abs(dt4_2i(:)))), colormap(gray(256))); axis xy;
% title('Phase-doubled then interpolated level 4 co-efficients: Magnitude');
% subplot(2,2,2); subimage(255*(mod(angle(dt{3}),2*pi)/(2*pi)), colormap(hsv(256))); axis xy;
% title('Phase-doubled then interpolated level 4 co-efficients: Phase');
% subplot(2,2,3); subimage(complex2rgb(dt4_2i)); axis xy;
% title('Phase-doubled  then interpolated level 4 co-efficients: Magnitude and phase');
% subplot(2,2,4); quiver(real(dt4_2i), imag(dt4_2i)); set(gca, 'ytick', 1:8, 'yticklabel', -16:8:40);
% title('Phase-doubled  then interpolated level 4 co-efficients: Magnitude and phase');

%
figure; 
subplot(2,2,1); subimage(255*(abs(ilp34_i2)/max(abs(ilp34_i2(:)))), colormap(gray(256))); axis xy;
title('Interpolated then phase-doubled inter-level co-efficients: Magnitude');
subplot(2,2,2); subimage(255*(mod(angle(ilp34_i2),2*pi)/(2*pi)), colormap(hsv(256))); axis xy;
title('Interpolated then phase-doubled inter-level co-efficients: Phase');
subplot(2,2,3); subimage(complex2rgb(ilp34_i2)); axis xy;
title('Interpolated  then phase-doubled inter-level co-efficients: Magnitude and phase');
subplot(2,2,4); quiver(real(ilp34_i2), imag(ilp34_i2));
title('Interpolated  then phase-doubled inter-level co-efficients: Magnitude and phase'); set(gca, 'ytick', 1:8, 'yticklabel', -16:8:40);

% figure; 
% subplot(2,2,1); subimage(255*(abs(ilp34_2i)/max(abs(ilp34_2i(:)))), colormap(gray(256))); axis xy;
% title('Phase-doubled then interpolated inter-level co-efficients: Magnitude');
% subplot(2,2,2); subimage(255*(mod(angle(ilp34_2i),2*pi)/(2*pi)), colormap(hsv(256))); axis xy;
% title('Phase-doubled then interpolated inter-level co-efficients: Phase');
% subplot(2,2,3); subimage(complex2rgb(ilp34_2i)); axis xy;
% title('Phase-doubled  then interpolated inter-level co-efficients: Magnitude and phase');
% subplot(2,2,4); quiver(real(ilp34_2i), imag(ilp34_2i));
% title('Phase-doubled  then interpolated inter-level co-efficients: Magnitude and phase'); set(gca, 'ytick', 1:8, 'yticklabel', -16:8:40);
%%

%--------------------------------------------------------------------------
% Other stuff that is prob not important now
%--------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I think the problem is where the interpolated data points are sampled. So
% lets test creating ilp coefficients for a set of offset sampled points
dim = size(dt{4}, 1);
for offset = -1:0.2:1
    
    dt4_i = zeros(size(dt{3}));

    %interpolate level 4 coeffs
    for ii = 1:16
        dt_mag = interp1(1:dim, abs(dt{4}(:,ii)), (0.5:0.5:dim)+offset, 'linear', 'extrap');
        dt_phase = interp1(1:dim, angle(dt{4}(:,ii)), (0.5:0.5:dim)+offset, 'linear', 'extrap');
        dt4_i(:,ii) = dt_mag .* exp(i*dt_phase);  %#ok       
    end
    
    %phase double them
    dt4_i2 =  abs(dt4_i).*complex(cos(2*angle(dt4_i)), sin(2*angle(dt4_i)));
    
    %compute and display the ilp coeffs
    ilp34_i2 = dt{3}.*conj(dt4_i2);
    figure; quiver(real(ilp34_i2), imag(ilp34_i2));
    title(['Inter-level product co-efficients: Interpolated then phase double. Offset = ', num2str(offset)] ); set(gca, 'ytick', 1:8, 'yticklabel', -16:8:40);
end
    
%%
%%%%%%%%%%%%%%%%%%%%%%
% Display results
%%%%%%%%%%%%%%%%%%%%%%

% Display step-edges
figure; imagesc(step_edges1D); colormap(gray); axis xy; set(gca, 'ytick', 8:8:64, 'yticklabel', -16:8:40);
hold on;
plot(repmat([0.5 16.5]', 1, 4), [(8:16:56); (8:16:56)], 'r:');
title('Step-edges: Level 4 sampling points marked');

%Display DT coefficients as phasors
figure; quiver(real(dt{3}), imag(dt{3})); set(gca, 'ytick', 1:8, 'yticklabel', -16:8:40);
title('Level 3 co-efficients');

figure; quiver(real(dt{4}), imag(dt{4}), 0); set(gca, 'ytick', 1:0.5:4.5, 'yticklabel', -16:8:40); axis equal
hold on; quiver(real(dt4_2), imag(dt4_2), 0, 'r'); set(gca, 'ytick', 1:0.5:4.5, 'yticklabel', -16:8:40); axis equal
title('Level 4 co-efficients: original in blue, phase doubled in red');

figure; quiver(real(dt4_i), imag(dt4_i)); set(gca, 'ytick', 1:8, 'yticklabel', -16:8:40);
title('Level 4 co-efficients: Interpolated');

figure; quiver(real(dt4_2i), imag(dt4_2i)); set(gca, 'ytick', 1:8, 'yticklabel', -16:8:40); 
title('Level 4 co-efficients: Phase doubled then interpolate');
figure; quiver(real(dt4_i2), imag(dt4_i2)); set(gca, 'ytick', 1:8, 'yticklabel', -16:8:40);
title('Level 4 co-efficients: Interpolated then phase double');

figure; quiver(real(ilp34_2i), imag(ilp34_2i)); set(gca, 'ytick', 1:8, 'yticklabel', -16:8:40);
title('Inter-level product co-efficients: Phase doubled then interpolated');
figure; quiver(real(ilp34_i2), imag(ilp34_i2));
title('Inter-level product co-efficients: Interpolated then phase double'); set(gca, 'ytick', 1:8, 'yticklabel', -16:8:40);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Maybe we should try working backwards instead
%Try subtracting pi/4 phase to each coefficient in level 3
dt3_pi4 = abs(dt{3}).*complex(cos(angle(dt{3}) - pi/4), sin(angle(dt{3}) - pi/4));
test_ilp = dt{3}.*conj(dt3_pi4);

figure; quiver(real(dt3_pi4), imag(dt3_pi4)); set(gca, 'ytick', 1:8, 'yticklabel', -16:8:40);
title('Level 3 co-efficients with pi/4 subtracted from phase');

figure; quiver(real(test_ilp), imag(test_ilp));
title('Inter-level product co-efficients: Well not really...'); set(gca, 'ytick', 1:8, 'yticklabel', -16:8:40);


%%
% Phase half the pi/4 subtracted level 3 coeffs - these should be
% equivalent to the interpolated level 4 coeffs
dt3_pi4_half = abs(dt3_pi4).*complex(cos(0.5*angle(dt3_pi4)), sin(0.5*angle(dt3_pi4)));
figure; quiver(real(dt3_pi4_half), imag(dt3_pi4_half));
title('Level co-efficients: Pi/4 subtracted then phase-halved'); set(gca, 'ytick', 1:8, 'yticklabel', -16:8:40);
%%
figure; 
subplot(2,1,1); imagesc(real(dt3_pi4_half));
subplot(2,1,2); imagesc(real(dt{4}));
figure; 
subplot(2,1,1); imagesc(imag(dt3_pi4_half));
subplot(2,1,2); imagesc(imag(dt{4}));

%%
% No do the opposite interpolation... that is for each of the real and
% imaginary parts, look at what location the non-interpolated level 4
% coefficients would have to have to generate the features in dt3_pi4_half
for ii = 1:16    
    real_pos(:,ii) = zeros(size(dt{4}, 1), 1); %#ok
    imag_pos(:,ii) = zeros(size(dt{4}, 1), 1); %#ok
    real_pos(:,ii) = interp1(real(dt3_pi4_half(:,ii)), 1:8, real(dt{4}(:,ii)), 'linear', 'extrap'); %#ok
    imag_pos(:,ii) = interp1(imag(dt3_pi4_half(:,ii)), 1:8, imag(dt{4}(:,ii)), 'linear', 'extrap'); %#ok
end


%%
%Look at phase of the test coefficients
figure;
for coeff = 1:8
    subplot(4, 2, coeff); hold on; xlabel(['Coeff: ', num2str(coeff)]);
    plot(1:16, angle(dt3_pi4_half(coeff,:)), 'b');
    plot(1:16, angle(dt3_pi4_half(coeff,:)), 'rx');
    plot([1 16], [0 0], 'r-');
end
%%
figure;
for coeff = 1:8
    subplot(4, 2, coeff); hold on; xlabel(['Coeff: ', num2str(coeff)]);
    plot(1:16, angle(dt4_i(coeff,:)), 'b');
    plot(1:16, angle(dt4_i(coeff,:)), 'rx');
    plot([1 16], [0 0], 'r-');
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A new interpolation attempt...
%Interpolate the 4th level co-efficients




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A new interpolation attempt...
%Interpolate the 4th level co-efficients
for ii = 1:16    
    dt4_j(:,ii) = zeros(size(dt{3}, 1), 1); %#ok
    real_part = interp1(1.5:2:7.5, real(dt{4}(:,ii)), 1:8, 'linear', 'extrap');
    imag_part = interp1(2:2:8, imag(dt{4}(:,ii)), 1:8, 'linear', 'extrap');
    
    dt4_j(:,ii) = complex(real_part, imag_part); %#ok
end
figure; quiver(real(dt4_j), imag(dt4_j));
title('Level 4 co-efficients: New method of interpolation'); set(gca, 'ytick', 1:8, 'yticklabel', -16:8:40);
%%
%Phase double the interpolated 4th level co-efficients
dt4_j2 =  abs(dt4_j).*complex(cos(2*angle(dt4_j)), sin(2*angle(dt4_j)));
ilp_j = dt{3}.*conj(dt4_j2);
figure; quiver(real(ilp_j), imag(ilp_j));
title('ILP co-efficients: New method of interpolation'); set(gca, 'ytick', 1:8, 'yticklabel', -16:8:40);