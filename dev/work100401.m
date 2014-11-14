x = 1:5;
z = [2 5 3 4 1];

A1 = trapz(z);
A2 = z(1)/2 + z(2) + z(3) + z(4) + z(5)/2;

xi = 1:0.1:5;
zi = interp1(x, z, xi, 'linear');

figure; hold on;
plot(x,z);
plot(xi, zi, 'rx');

A3 = conv2(z, [.5 1 1 1 .5], 'valid');

A4 = trapz(zi)*.1;

B1 = trapz(zi(6:end)) * .1;
B2 = (0.5*(zi(6) + z(2))/ 2) + z(2)/2 + z(3) + z(4) + z(5)/2;
B3 = conv2(z, fliplr([.5*(x(2)-xi(6)).^2 .5*((x(2)-xi(6))*(xi(6)-x(1)+ 1) + 1) 1 1 .5]), 'valid');
B4 = .5*(x(2)-xi(6)).^2 *z(1) + .5*(2*(x(2)-xi(6)) - (x(2)-xi(6)).^2 + 1)*z(2) + z(3) + z(4) + z(5)/2;

%%
clc
for i = 1:9
    (.5*(x(2)-xi(i))*(zi(i) + z(2)) + .5*z(2) ) - ...
        ( .5*(x(2)-xi(i)).^2 *z(1) + .5*(2*(x(2)-xi(i)) - (x(2)-xi(i)).^2 + 1)*z(2))

%         .5*z(1)*(x(2)-xi(i)).^2 + .5*z(2)*(xi(i)-x(1))*(x(2)-xi(i)) + .5*z(2)
%     .5*(x(2)-xi(i))*(z(1) + (z(2)-z(1))*(xi(i)-x(1)) + z(2))
%     .5*(x(2)-xi(i))*(z(1)*(1 - xi(i) + x(1)) + z(2)*(xi(i)-x(1)) + z(2))
%     .5*(x(2)-xi(i))*(z(1)*(x(2) - xi(i)) + z(2)*(xi(i)-x(1)) + z(2))
    
end
%%
f_length = 11;
 
x1 = floor(xi);
x2 = ceil(xi);
c = (f_length - 1) / 2;

f(1:x1-1) = 0;
f(x1) = (x2-xi).^2;
f(x2) = 2*(x2-x1) - (x2-xi).^2 + 1;
f(x2+1:c) = 2;

f = [f 2 fliplr(f)];
%%
figure; hold all;
x = linspace(0, pi, 100);
for a = 0:0.1:1
    base_wave = sin(x); for ii = 3:2:11; base_wave = base_wave + a*sin(ii*x)/ii; end
    base_wave = base_wave / max(base_wave(:));
    %base_wave = 4*base_wave / (pi*a);
    plot(x, base_wave);
end