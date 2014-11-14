function Y = pi_rotate(X)
%Equivalent to Y = mb_mod(X, 2*pi);
%Transform the elements of X into the range [-pi, pi)
Y = mod(X, 2*pi);
Y(Y >= pi) = Y(Y >= pi) - 2*pi;