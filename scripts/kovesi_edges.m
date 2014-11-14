[x phi] = meshgrid(linspace(0,3*pi,181), linspace(0,pi,361));

waves = zeros(361,181);
p = 1;
figure; 
a1 = subplot(4,2,2); hold all; ylabel('0');
a2 = subplot(4,2,4); hold all; ylabel('pi/6');
a3 = subplot(4,2,6); hold all; ylabel('pi/3');
a4 = subplot(4,2,8); hold all; ylabel('pi/2');
for n = 0:181
    
    wave = sin((2*n + 1)*x + phi) / (2*n + 1)^p;
    waves = waves + wave;
    if n < 5
        plot(a1, x(1,:), wave(1,:), '--');
        plot(a2, x(61,:), wave(61,:), '--');
        plot(a3, x(121,:), wave(121,:), '--');
        plot(a4, x(181,:), wave(181,:), '--');
    end
end
plot(a1, x(1,:), waves(1,:), 'k');
plot(a2, x(61,:), waves(61,:), 'k');
plot(a3, x(121,:), waves(121,:), 'k');
plot(a4, x(181,:), waves(181,:), 'k');

subplot(4,2, [1 3 5 7]); imagesc(waves); axis image; colormap(gray(256)); hold on;
title('Interpolation of a step feature to a line feature');
set(gca,...
    'xtick', [],...
    'xticklabel', [],...
    'ytick', [1 61 121 181],...
    'yticklabel', {'0', 'pi/6', 'pi/3', 'pi/2'});
plot([0 0 0 0; 1 1 1 1]*361, [1 61 121 181; 1 61 121 181], 'r--');

[g, dg, ddg] = gaussian_filters_1d(2);
%G2 = ddg;
%H2 = imag(hilbert(G2));
%H2 = H2 - mean(H2(:));
G2 = g' * ddg;
%H2 = g' * dg;
H2 = G2;
for ii = 1:size(H2,1)
    H2(ii,:) = imag(hilbert(H2(ii,:)));
end
I_g2 = imfilter(waves, G2, 'replicate');
I_h2 = imfilter(waves, H2, 'replicate');

figure;
subplot(1,2,1);imagesc(I_g2); axis image; colormap(gray(256));
title('Response to G_2');
subplot(1,2,2);imagesc(I_h2); axis image; colormap(gray(256));
title('Response to H_2');
figure; hold on;
plot(I_g2(:,61), I_h2(:,61), 'b-'); axis equal;
plot(I_g2(:,121), I_h2(:,121), 'r-'); axis equal;
xlabel('G_2');
ylabel('H_2');
title({'Filter response to G_2 and H_2 sampled vertically along features'; 'Blue response from left feature, red responses from right feature'});

figure; hold on;
plot([0 pi; 0 pi; 0 pi; 0 pi; 0 pi]', [pi pi; pi/2 pi/2; 0 0; -pi/2 -pi/2; -pi -pi]', 'k--');
plot(linspace(0,pi,361), atan2(I_h2(:,61), I_g2(:,61)), 'b-');
plot(linspace(0,pi,361), atan2(I_h2(:,121), I_g2(:,121)), 'r-');
set(gca,...
    'ytick', [-pi -pi/2 0 pi/2 pi],...
    'yticklabel', {'-pi', '-pi/2', '0', 'pi/2', 'pi'});
title({'Phase of response (i.e. arctan(H_2/G_2) ) sampled vertically along features'; 'Blue response from left feature, red responses from right feature'});

figure; hold on;
plot(linspace(0,pi,361), sqrt(I_h2(:,61).^2 + I_g2(:,61).^2), 'b-');
plot(linspace(0,pi,361), sqrt(I_h2(:,121).^2 + I_g2(:,121).^2), 'r--');
title({'Magnitude of response (i.e. sqrt(H_2^2 + G_2^2) ) sampled vertically along features'; 'Blue response from left feature, red responses from right feature'});
set(gca, 'ylim', [0 max(get(gca, 'ylim'))]);