% Theoretical work on the DT-CWT

%Generate some step_edges
base_edge(:,:,1) = zeros(64);
base_edge(31,:) = 63;
base_edge(32,:) = 127;
base_edge(33,:) = 191;
base_edge(34:64,:) = 255;

test_edge = zeros(32,32,10);

for ii = 1:10
    temp = imrotate(base_edge, ii*18, 'crop');
    test_edge(:,:,ii) = temp(17:48, 17:48);
    clear temp;
    figure; image(test_edge(:,:,ii)); colormap(gray(256)); axis image;
end

mammogram = u_load('C:\isbe\mammograms\new_CAD\bMP_2004_half\o04_001LCC.mat');
mammo_patch = double(mammogram(1000:1031, 1000:1031)); clear mammogram;
mammo_edge = zeros(32,32,10);

for ii = 1:10
    mammo_edge(:,:,ii) = mammo_patch + test_edge(:,:,ii)/10;
    figure; image(mammo_edge(:,:,ii)); colormap(gray(256)); axis image;
end
%%
%Generate some line_edges
base_line = zeros(32);
base_line(:,15) =40;
base_line(:,18) =40;
base_line(:,16:17) =60;

mammo_line = zeros(32,32,10);

for ii = 1:10;
    mammo_line(:,:,ii) = mammo_patch + imrotate(base_line, ii*18, 'crop');
    figure; imagesc(mammo_line(:,:,ii)); axis image; colormap(gray(256));
end

mammo_star = mammo_patch;
base_line2 = base_line;
base_line2(12:20,:) = 0;

for ii = 1:4;
    mammo_star = mammo_star + imrotate(base_line2, ii*45, 'crop');
end
figure; imagesc(mammo_star); axis image; colormap(gray(256));

%Play around with some dual_trees
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step edges
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[dt_edge sc_edge] = dtwavexfm2(mammo_edge(:,:,1), 3, 'near_sym_b','qshift_b');
axis_lims = [1 16 1 16 min(min(real(dt_edge{1}(:))), min(imag(dt_edge{1}(:)))) max(max(real(dt_edge{1}(:))), max(imag(dt_edge{1}(:))))];
for ori = 1:6
    figure; 
    h1 = subplot(1,2,1); mesh(real(dt_edge{1}(:,:,ori))); axis(axis_lims);
    h2 = subplot(1,2,2); mesh(imag(dt_edge{1}(:,:,ori))); axis(axis_lims);
    hlinks(ori) = linkprop([h1 h2], {'CameraPosition','CameraUpVector'}); %#ok
end
%
axis_lims = [1 16 1 16 min(abs(dt_edge{1}(:))) max(abs(dt_edge{1}(:)))];

figure;
for ori = 1:6
     
    ho(ori) = subplot(2,3,ori); mesh(abs(dt_edge{1}(:,:,ori))); axis(axis_lims); %#ok
    
end
ho_links = linkprop(ho, {'CameraPosition','CameraUpVector'});
%%
figure; quiver(real(dt_edge{1}(:,:,4)), imag(dt_edge{1}(:,:,4))); axis image;
figure; quiver(real(dt_edge{2}(:,:,4)), imag(dt_edge{2}(:,:,4))); axis image;
%%
[xx yy] = meshgrid(1:16,1:16);
dt_2_4_int = interp2(1:2:15,1:2:15,dt_edge{2}(:,:,4),xx,yy, '*linear', 0);
%%
figure;
h1 = subplot(1,2,1); mesh(abs(dt_edge{2}(:,:,4)));
h2 = subplot(1,2,2); mesh(abs(dt_2_4_int));
hlinks(ori) = linkprop([h1 h2], {'CameraPosition','CameraUpVector'}); %#ok
%%
phase_diff = dt_edge{1}(:,:,4) - dt_2_4_int;
figure;
quiver(real(phase_diff), imag(phase_diff)); axis image;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[dt_star sc_star] = dtwavexfm2(mammo_star, 3, 'near_sym_b','qshift_b');
%%
axis_lims = [1 16 1 16 min(min(real(dt_star{1}(:))), min(imag(dt_star{1}(:)))) max(max(real(dt_star{1}(:))), max(imag(dt_star{1}(:))))];
for ori = 1:6
    figure; 
    h1 = subplot(1,2,1); mesh(real(dt_star{1}(:,:,ori))); axis(axis_lims);
    h2 = subplot(1,2,2); mesh(imag(dt_star{1}(:,:,ori))); axis(axis_lims);
    hlinks(ori) = linkprop([h1 h2], {'CameraPosition','CameraUpVector'}); %#ok
end

%%
axis_lims = [1 16 1 16 0 max(abs(dt_star{1}(:)))];
figure;
for ori = 1:6
    h(ori) = subplot(2,3,ori); mesh(abs(dt_star{1}(:,:,ori))); axis(axis_lims); %#ok 
    
end
hlinks_ori(end+1) = linkprop(h, {'CameraPosition','CameraUpVector'});

%%
mammogram = u_load('C:\isbe\mammograms\new_CAD\bMP_2004_half\o04_001LCC.mat');
mammo_patch64 = double(mammogram(1000:1063, 1000:1063)); clear mammogram;
%%
base_line3 = zeros(64);
base_line3(1:24,31) =0;
base_line3(1:24,34) =0;
base_line3(1:24,32:32) =60;

mammo_star64 = mammo_patch64;

angles = [75 45 15 -15 -45 -70];
%%
for ii = 1:6;
    mammo_star64 = mammo_star64 + imrotate(base_line3, angles(ii), 'crop') + imrotate(base_line3, angles(ii)+180, 'crop');
end
figure; imagesc(mammo_star64); axis image; colormap(gray(256));
%%
[dt_star64 sc_star64] = dtwavexfm2(mammo_star64, 3, 'near_sym_b','qshift_b');
%%
axis_lims = [1 32 1 32 0 max(abs(dt_star64{1}(:)))];
figure;
for ori = 1:6
    h(ori) = subplot(2,3,ori); mesh(abs(dt_star64{1}(:,:,ori))); axis(axis_lims); %#ok 
    
end
hlinks_ori(end+1) = linkprop(h, {'CameraPosition','CameraUpVector'});
%%
axis_lims = [1 32 1 32 min(real(dt_star64{1}(:))) max(real(dt_star64{1}(:)))];
figure;
for ori = 1:6
    h(ori) = subplot(2,3,ori); mesh(real(dt_star64{1}(:,:,ori))); axis(axis_lims); %#ok 
    
end
hlinks_ori(end+1) = linkprop(h, {'CameraPosition','CameraUpVector'});
%
axis_lims = [1 32 1 32 min(imag(dt_star64{1}(:))) max(imag(dt_star64{1}(:)))];
figure;
for ori = 1:6
    h(ori) = subplot(2,3,ori); mesh(imag(dt_star64{1}(:,:,ori))); axis(axis_lims); %#ok 
    
end
hlinks_ori(end+1) = linkprop(h, {'CameraPosition','CameraUpVector'});
%%
axis_lims = [1 32 1 32 0 max(abs(dt_star64{1}(:)))];
figure;
for ori = 1:6
    temp_real = real(dt_star64{1}(:,:,ori));
    temp_imag = imag(dt_star64{1}(:,:,ori));
    %temp_real(temp_real < 0) = 0;
    %temp_imag(temp_imag < 0) = 0;
    
    h(ori) = subplot(2,3,ori); mesh(sqrt(temp_real.^2 + temp_imag.^2)); axis(axis_lims); %#ok 
    
end
hlinks_ori(end+1) = linkprop(h, {'CameraPosition','CameraUpVector'});
%%

figure;
for ori = 1:6
    subplot(2,3,ori); imagesc(abs(dt_star64{1}(:,:,ori))); axis image; %#ok 
    
end
%%
axis_lims = [1 16 1 16 0 max(abs(dt_star64{2}(:)))];
figure;
for ori = 1:6
    h(ori) = subplot(2,3,ori); mesh(abs(dt_star64{2}(:,:,ori))); axis(axis_lims); %#ok 
    
end
hlinks_ori(end+1) = linkprop(h, {'CameraPosition','CameraUpVector'});
%
figure;
for ori = 1:6
    subplot(2,3,ori); imagesc(abs(dt_star64{2}(:,:,ori))); axis image; %#ok 
    
end
%%
axis_lims = [1 16 1 16 min(min(real(dt_star{2}(:))), min(imag(dt_star{2}(:)))) max(max(real(dt_star{2}(:))), max(imag(dt_star{2}(:))))];
figure;
for ori = 1:6
    h(ori) = subplot(2,3,ori);  %#ok
    mesh(real(dt_star64{2}(:,:,ori)), 'FaceColor', 'none', 'EdgeColor', 'b'); hold on;
    mesh(imag(dt_star64{2}(:,:,ori)), 'FaceColor', 'none', 'EdgeColor', 'r'); axis(axis_lims); 
    
end
hlinks_ori(end+1) = linkprop(h, {'CameraPosition','CameraUpVector'});
%%
axis_lims = [1 32 1 32 min(min(real(dt_star64{1}(:))), min(imag(dt_star64{1}(:)))) max(max(real(dt_star64{1}(:))), max(imag(dt_star64{1}(:))))];
figure;
for ori = 1:6
    h(ori) = subplot(2,3,ori);  %#ok
    mesh(real(dt_star64{1}(:,:,ori)), 'FaceColor', 'none', 'EdgeColor', 'b'); hold on;
    mesh(imag(dt_star64{1}(:,:,ori)), 'FaceColor', 'none', 'EdgeColor', 'r'); axis(axis_lims); 
    
end
hlinks_ori(end+1) = linkprop(h, {'CameraPosition','CameraUpVector'});
%%
figure;
for ori = 1:6
    subplot(2,3,ori); imagesc(abs(dt_star64{2}(:,:,ori))); axis image; %#ok 
    
end
%%
figure;
for ori = 1:6
    subplot(2,3,ori); hist(abs(reshape(dt_star64{1}(:,:,ori), [],1)), 30); %#ok 
    
end
%%
figure;
for ori = 1:6
    subplot(2,3,ori); imagesc(abs(dt_star64{1}(:,:,ori))>10); axis image; 
    
end
%%
for ori = 1:6
    figure('name', ['Orientation = ', num2str(angles(ori))]);
    subplot(1,2,1); imagesc(real(dt_star64{1}(:,:,ori))); axis image;
    subplot(1,2,2); imagesc(imag(dt_star64{1}(:,:,ori))); axis image;    
end

%%
[dt_edge2 sc_edge2] = dtwavexfm2(mammo_edge(:,:,2), 3, 'near_sym_b','qshift_b');
figure;
for ori = 1:6
    subplot(2,3,ori); imagesc(abs(dt_edge2{1}(:,:,ori))); axis image; %#ok 
    
end
for ori = 1:6
    subplot(2,3,ori); imagesc(abs(dt_edge2{2}(:,:,ori))); axis image; %#ok 
    
end