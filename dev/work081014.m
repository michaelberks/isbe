% Work on Peter Kovesi's paper: Edges are not just steps

%Create three Gating images, figure 2; from the paper


x = linspace(0, 3*pi, 256);

gatings = zeros(256,256,3);
p = [0.5, 1.0, 1.5];
pheta = linspace(0, pi/2, 256);
lim = 256;
%%
for i = 1:length(pheta);
    
    for g = 1:3;
        gatings(i,:,g) = fourier_sum(x, p(g), pheta(i), lim, 0);
    end    
end
%
c_lims(1) = min(gatings(:));
c_lims(2) = max(gatings(:));
%
for g = 1:3
    figure; 
    subplot(4,2,[1 3 5 7]); imagesc(gatings(:,:,g)); axis image; colormap(gray); caxis(c_lims);
    subplot(4,2,2); plot(x, gatings(1,:,g), 'r');
    subplot(4,2,4); plot(x, gatings(86,:,g), 'g');
    subplot(4,2,6); plot(x, gatings(171,:,g), 'b');
    subplot(4,2,8); plot(x, gatings(256,:,g), 'm');
end
%%

%working with the 2nd gating, convert to animage in range 0:256

gate_image = gatings(:,:,2);
gate_image = gate_image - min(gate_image(:));
gate_image = gate_image / max(gate_image(:));
gate_image = gate_image*255;

%%
[dt_gate sc_gate] = dtwavexfm2(gate_image, 3, 'near_sym_b','qshift_b');

axis_lims = [1 128 1 128 min(min(real(dt_gate{1}(:))), min(imag(dt_gate{1}(:)))) max(max(real(dt_gate{1}(:))), max(imag(dt_gate{1}(:))))];
for ori = 1:6
    figure; 
    h1 = subplot(1,2,1); mesh(real(dt_gate{1}(:,:,ori))); axis(axis_lims);
    h2 = subplot(1,2,2); mesh(imag(dt_gate{1}(:,:,ori))); axis(axis_lims);
    hlinks(ori) = linkprop([h1 h2], {'CameraPosition','CameraUpVector'}); %#ok
end
%%
axis_lims = [1 128 1 128 min(abs(dt_gate{1}(:))) max(abs(dt_gate{1}(:)))];

figure;
for ori = 1:6
     
    ho(ori) = subplot(2,3,ori); mesh(abs(dt_gate{1}(:,:,ori))); axis(axis_lims); %#ok
    
end
ho_links = linkprop(ho, {'CameraPosition','CameraUpVector'});
%%
figure;
for ori = 3:4
     
    ao(ori) = subplot(1,2,ori-2); quiver(real(dt_gate{1}(1:4:end,1:4:end,ori)), real(dt_gate{1}(1:4:end,1:4:end,ori))); %#ok
    
end
linkaxes(ao);