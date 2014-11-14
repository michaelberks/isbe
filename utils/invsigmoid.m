function x = invsigmoid(y)

% Avoid infinities by shrinking toward 0.5;
y = ((y - 0.5) * 0.999) + 0.5;
x = -log((1 - y) ./ y);
