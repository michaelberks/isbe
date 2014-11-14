function E_c = ft_2d_gradient(phi_x)

alpha = 0.0005;
lamda = 0.5;


[gx gy] = gradient(phi_x);

gU_phi_x = lamda * (phi_x.^3 - phi_x) + alpha * (1 - phi_x.^2);

E_O_phi_x = 0.5*(gx.^2 + gy.^2) + gU_phi_x;
E_O = sum(E_O_phi_x(:));

E_S = CalculateEquivalentPhaseFieldEnergy(gx, gy);

E_c = E_O + E_S;

