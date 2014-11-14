figure; hold on;

mesh(2:16, 5*(1:36), oe_old, 'edgecolor', 'r', 'facecolor', 'r', 'facealpha', 0.5); %shading interp;
mesh(2:2:16, 5*(1:36), oe_new, 'edgecolor', 'g', 'facecolor', 'g', 'facealpha', 0.5); %shading interp;
mesh(2:2:16, 5*(1:36), oe_g2, 'edgecolor', 'b', 'facecolor', 'b', 'facealpha', 0.5); %shading interp;
%
figure; hold on;
% plot(repmat((2:2:16)',1,36), oe_new', 'g:');
% plot(repmat((2:2:16)',1,36), oe_g2', 'b:');
% plot(repmat((2:16)',1,36), oe_old', 'r:');

plot(2:2:16, mean(oe_new,1), 'k', 'linewidth', 2);
plot(2:2:16, mean(oe_g2,1), 'k', 'linewidth', 2);
plot(2:16, mean(oe_old,1), 'k', 'linewidth', 2);
plot(2:2:16, mean(oe_new,1), 'g--', 'linewidth', 2);
plot(2:2:16, mean(oe_g2,1), 'b--', 'linewidth', 2);
plot(2:16, mean(oe_old,1), 'r--', 'linewidth', 2);
%
figure; hold on;
plot(repmat((1:36)',1,8), oe_new, 'g:');
plot(repmat((1:36)',1,8), oe_g2, 'b:');
plot(repmat((1:36)',1,15), oe_old, 'r:');

plot(1:36, mean(oe_new,2), 'k', 'linewidth', 4);
plot(1:36, mean(oe_g2,2), 'k', 'linewidth', 4);
plot(1:36, mean(oe_old,2), 'k', 'linewidth', 4);
plot(1:36, mean(oe_new,2), 'g--', 'linewidth', 4);
plot(1:36, mean(oe_g2,2), 'b--', 'linewidth', 4);
plot(1:36, mean(oe_old,2), 'r--', 'linewidth', 4);
