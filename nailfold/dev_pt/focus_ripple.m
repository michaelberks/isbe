clc; 

n = 301;
k = 20;
gt = zeros(3*k, n);
gt_centre = 2*k - k*sin(t);

t = linspace(0,pi,n);
for i = 1:n
    gt(floor(gt_centre(i)), i) = 1;
end

f = isbe_fspecial('gaussian', [15,15], 4);
gt = conv2(gt, f, 'same');
gt = gt / max(gt(:));

gt = gt + 0.2*randn(size(gt));

% Initial condition
m = 15;

vel = cos( linspace(0, 2*pi, m+1) );
vel = vel(1:end-1);

z = zeros(1, m);

% Do initial scan without moving
x = 1;
y0 = 2*k;
mi = 1;
for i = 1:m
    z(mi) = interp2(gt, x, y0, '*linear');
    y0 = y0 + vel(mi);
    mi = rem(mi, m) + 1;
end
bias = vel(:)'*z(:) / sum(z(:));

figure(1); clf; hold on; colormap(gray(256));
    imagesc(gt);
    axis('image', 'ij');
    plot(gt_centre, 'b.');

mi = 1;
x = linspace(1, n, 2*n);
y = size(x);
y(1) = y0;
for i = 2:length(x)
    y(i) = y(i-1) + 2*vel(mi) + 2*bias;
    z(mi) = interp2(gt, x(i), y(i), '*linear');
    zn = z(:).^2;
    bias = vel(:)' * ( zn(:) / sum(zn(:)) );
        
    mi = rem(mi, m) + 1;
end

figure(1);
    plot(x, y, 'r.');