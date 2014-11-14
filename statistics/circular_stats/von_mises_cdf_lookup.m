function [P_lookup] = von_mises_cdf_lookup()

px_per_mm = 100/9;
r_max = 400;
r_min = floor(px_per_mm*5);
R = px_per_mm*4;

mu_i = linspace(0, pi, 180);
R_j = linspace(0, 1, 100);
r_k = r_min:r_max;

n_i = length(mu_i);
n_j = length(R_j);
n_k = length(r_k);
P_lookup = zeros(n_i, n_j, n_k);

for ii = 1:n_i
    for jj = 1:n_j
        for kk = 1:n_k
            
            alpha = pi * R / r_k(kk);
            mu = mu_i(ii);
            kappa = Ainv(R_j(jj)); 
            P_lookup(ii, jj, kk) = von_mises_cdf([-alpha alpha], mu, kappa);
        end
    end
    display(['Complete for ii = ' num2str(ii)]);
end
save([asymmetryroot 'data/misc/von_mises_P_lookup.mat'], 'P_lookup');

function k = Ainv(R)
% Compute the inverse of the function A(k) = I1(k)/I0(k) where Ij are
% modified bessel functions of the first kind

%Choose a set of k values to act as interpolation points
kk = linspace(0, 100, 1e5);
I_0 = besseli(0, kk);
I_1 = besseli(1, kk);

A_k = I_1 ./ I_0;

k = interp1(A_k, kk, R, 'linear', 'extrap');
