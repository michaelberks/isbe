function m = mb_mod(x, y)
% I want a function that behaves like mod, but returns me m = mod(x,y) such
% that -y < m < y and abs(m) in minimal. For example, using matlab's mod
% function mod(5,3) = 2, but I want to return -1

y_half = abs(y / 2);

m = mod(x, y);

m(abs(m) > y_half) = m(abs(m) > y_half) - y;
