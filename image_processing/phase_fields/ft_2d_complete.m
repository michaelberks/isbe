function E_c = ft_2d_complete(phi_x, d, alpha, lambda)

if nargin < 2
    d = 3;
    alpha = 0.0005;
    lambda = 0.5;
end

[gx gy] = gradient(phi_x);

E_O = CalculateEO(phi_x, gx, gy, alpha, lambda);

E_S = CalculateEquivalentPhaseFieldEnergy(gx, gy, d);

E_c = E_O + E_S;

