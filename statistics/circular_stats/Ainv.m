function k = Ainv(R)
% Compute the inverse of the function A(k) = I1(k)/I0(k) where Ij are
% modified bessel functions of the first kind

%Choose a set of k values to act as interpolation points
kk = linspace(0, 100, 1e5);
I_0 = besseli(0, kk);
I_1 = besseli(1, kk);

A_k = I_1 ./ I_0;

k = interp1(A_k, kk, R, 'linear', 'extrap');