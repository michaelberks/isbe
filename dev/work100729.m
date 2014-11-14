y = complex(randn(4,1), randn(4,1));
x = 0:3;
xi = 0.01:0.01:2.99;
%--------------------------------------------------------------------------
%%
method = 'linear';

yi1 = interp1(x, y, xi, method);
yi2 = complex(interp1(x, real(y), xi, method), interp1(x, imag(y), xi, method));
yi3 = interp1(x, abs(y), xi, method) .* exp(i*interp1(x, angle(y), xi, method));
yi4 = interp1(x, abs(y), xi, method) .* exp(i*angle(interp1(x, y ./ (abs(y)+1e-6), xi, 'linear')));
%
figure;
subplot(3,1,1); hold on;
plot3([-1 4], [0 0], [0 0], 'k');
for ii = 1:length(xi)
    plot3([xi(ii) xi(ii)], [0 real(yi1(ii))], [0 imag(yi1(ii))], 'r', 'linewidth', 2);
    plot3(xi(ii), real(yi1(ii)), imag(yi1(ii)), 'ro');
end
for ii = 1:length(x)
    plot3([x(ii) x(ii)], [0 real(y(ii))], [0 imag(y(ii))], 'b', 'linewidth', 2);
    plot3(x(ii), real(y(ii)), imag(y(ii)), 'bo');
end
subplot(3,1,2); hold on;
plot3([-1 4], [0 0], [0 0], 'k');
for ii = 1:length(xi)
    plot3([xi(ii) xi(ii)], [0 real(yi2(ii))], [0 imag(yi2(ii))], 'g', 'linewidth', 2);
    plot3(xi(ii), real(yi2(ii)), imag(yi2(ii)), 'go');
end
for ii = 1:length(x)
    plot3([x(ii) x(ii)], [0 real(y(ii))], [0 imag(y(ii))], 'b', 'linewidth', 2);
    plot3(x(ii), real(y(ii)), imag(y(ii)), 'bo');
end
subplot(3,1,3); hold on;
plot3([-1 4], [0 0], [0 0], 'k');
for ii = 1:length(xi)
    plot3([xi(ii) xi(ii)], [0 real(yi3(ii))], [0 imag(yi3(ii))], 'm', 'linewidth', 2);
    plot3(xi(ii), real(yi3(ii)), imag(yi3(ii)), 'mo');
end
for ii = 1:length(x)
    plot3([x(ii) x(ii)], [0 real(y(ii))], [0 imag(y(ii))], 'b', 'linewidth', 2);
    plot3(x(ii), real(y(ii)), imag(y(ii)), 'bo');
end
%
figure; hold on;
plot3([-1 4], [0 0], [0 0], 'k');
for ii = 1:length(xi)
    plot3([xi(ii) xi(ii)], [0 real(yi1(ii))], [0 imag(yi1(ii))], 'r', 'linewidth', 2);
    plot3(xi(ii), real(yi1(ii)), imag(yi1(ii)), 'ro');
    
    plot3([xi(ii) xi(ii)], [0 real(yi2(ii))], [0 imag(yi2(ii))], 'g', 'linewidth', 2);
    plot3(xi(ii), real(yi2(ii)), imag(yi2(ii)), 'go');
    
    plot3([xi(ii) xi(ii)], [0 real(yi3(ii))], [0 imag(yi3(ii))], 'm', 'linewidth', 2);
    plot3(xi(ii), real(yi3(ii)), imag(yi3(ii)), 'mo');
end
for ii = 1:length(x)
    plot3([x(ii) x(ii)], [0 real(y(ii))], [0 imag(y(ii))], 'b', 'linewidth', 2);
    plot3(x(ii), real(y(ii)), imag(y(ii)), 'bo');
end
%
figure;
subplot(2,1,1); hold on;
plot(x, abs(y), 'bx');
plot(xi, abs(yi1), 'r');
plot(xi, abs(yi2), 'g');
plot(xi, abs(yi3), 'm');
plot(xi, abs(yi4), 'c');
subplot(2,1,2); hold on;
plot(x, angle(y), 'bx');
plot(xi, angle(yi1), 'r');
plot(xi, angle(yi2), 'g');
plot(xi, angle(yi3), 'm');
plot(xi, angle(yi4), 'c');
%--------------------------------------------------------------------------

method = 'cubic';

yi1 = interp1(x, y, xi, method);
yi2 = complex(interp1(x, real(y), xi, method), interp1(x, imag(y), xi, method));
yi3 = interp1(x, abs(y), xi, method) .* exp(i*interp1(x, angle(y), xi, method));
yi4 = interp1(x, abs(y), xi, method) .* exp(i*angle(interp1(x, y ./ (abs(y)+1e-6), xi, 'linear')));

figure;
subplot(3,1,1); hold on;
plot3([-1 4], [0 0], [0 0], 'k');
for ii = 1:length(xi)
    plot3([xi(ii) xi(ii)], [0 real(yi1(ii))], [0 imag(yi1(ii))], 'r', 'linewidth', 2);
    plot3(xi(ii), real(yi1(ii)), imag(yi1(ii)), 'ro');
end
for ii = 1:length(x)
    plot3([x(ii) x(ii)], [0 real(y(ii))], [0 imag(y(ii))], 'b', 'linewidth', 2);
    plot3(x(ii), real(y(ii)), imag(y(ii)), 'bo');
end
subplot(3,1,2); hold on;
plot3([-1 4], [0 0], [0 0], 'k');
for ii = 1:length(xi)
    plot3([xi(ii) xi(ii)], [0 real(yi2(ii))], [0 imag(yi2(ii))], 'g', 'linewidth', 2);
    plot3(xi(ii), real(yi2(ii)), imag(yi2(ii)), 'go');
end
for ii = 1:length(x)
    plot3([x(ii) x(ii)], [0 real(y(ii))], [0 imag(y(ii))], 'b', 'linewidth', 2);
    plot3(x(ii), real(y(ii)), imag(y(ii)), 'bo');
end
subplot(3,1,3); hold on;
plot3([-1 4], [0 0], [0 0], 'k');
for ii = 1:length(xi)
    plot3([xi(ii) xi(ii)], [0 real(yi3(ii))], [0 imag(yi3(ii))], 'm', 'linewidth', 2);
    plot3(xi(ii), real(yi3(ii)), imag(yi3(ii)), 'mo');
end
for ii = 1:length(x)
    plot3([x(ii) x(ii)], [0 real(y(ii))], [0 imag(y(ii))], 'b', 'linewidth', 2);
    plot3(x(ii), real(y(ii)), imag(y(ii)), 'bo');
end

%
figure; hold on;
plot3([-1 4], [0 0], [0 0], 'k');
for ii = 1:length(xi)
    plot3([xi(ii) xi(ii)], [0 real(yi1(ii))], [0 imag(yi1(ii))], 'r', 'linewidth', 2);
    plot3(xi(ii), real(yi1(ii)), imag(yi1(ii)), 'ro');
    
    plot3([xi(ii) xi(ii)], [0 real(yi2(ii))], [0 imag(yi2(ii))], 'g', 'linewidth', 2);
    plot3(xi(ii), real(yi2(ii)), imag(yi2(ii)), 'go');
    
%     plot3([xi(ii) xi(ii)], [0 real(yi3(ii))], [0 imag(yi3(ii))], 'm', 'linewidth', 2);
%     plot3(xi(ii), real(yi3(ii)), imag(yi3(ii)), 'mo');
end
for ii = 1:length(x)
    plot3([x(ii) x(ii)], [0 real(y(ii))], [0 imag(y(ii))], 'b', 'linewidth', 2);
    plot3(x(ii), real(y(ii)), imag(y(ii)), 'bo');
end

figure;
subplot(2,1,1); hold on;
plot(x, abs(y), 'bx');
plot(xi, abs(yi1), 'r');
plot(xi, abs(yi2), 'g');
plot(xi, abs(yi3), 'm');
plot(xi, abs(yi4), 'c');
subplot(2,1,2); hold on;
plot(x, angle(y), 'bx');
plot(xi, angle(yi1), 'r');
plot(xi, angle(yi2), 'g');
plot(xi, angle(yi3), 'm');
plot(xi, angle(yi4), 'c');
%%
figure; hold on;
plot3([-1 4], [0 0], [0 0], 'k');
for ii = 1:length(xi)
%     plot3([xi(ii) xi(ii)], [0 real(yi1(ii))], [0 imag(yi1(ii))], 'r', 'linewidth', 2);
%     plot3(xi(ii), real(yi1(ii)), imag(yi1(ii)), 'ro');
    
    plot3([xi(ii) xi(ii)], [0 real(yi2(ii))], [0 imag(yi2(ii))], 'g', 'linewidth', 2);
    plot3(xi(ii), real(yi2(ii)), imag(yi2(ii)), 'go');
    
    plot3([xi(ii) xi(ii)], [0 real(yi3(ii))], [0 imag(yi3(ii))], 'm', 'linewidth', 2);
    plot3(xi(ii), real(yi3(ii)), imag(yi3(ii)), 'mo');
end
for ii = 1:length(x)
    plot3([x(ii) x(ii)], [0 real(y(ii))], [0 imag(y(ii))], 'b', 'linewidth', 2);
    plot3(x(ii), real(y(ii)), imag(y(ii)), 'bo');
end
%%
mammo = double(imread('C:\isbe\asymmetry_project\data\misc\o04_010RCC_1024_3797_3365.bmp'));
dt = dtwavexfm2(mammo, 4);
full_dt = dt_to_full_image(dt, 4);
%%
