function scale_model = generate_scale_model(mass_model)
% function to take an existing mass model and modify so the combined shape
% and texture parameters are conditioned on the scale parameter
scale_model = mass_model;

k_shape = length(scale_model.L_shape);
k_tex   = length(scale_model.L_tex);

combined_data = [scale_model.W_shape*scale_model.B_shape; scale_model.W_tex*scale_model.B_tex]';

scale_model.combined_data = combined_data;

[mean_com, P_com, B_com, L_com] = pca(combined_data, 0.98);

scale_model.mean_com = mean_com;
scale_model.P_com = P_com;
scale_model.B_com = B_com;
scale_model.L_com = L_com;

k_com = length(scale_model.L_com);

C_com = cov([B_com; scale_model.W_scale*scale_model.B_scale]');
C_com_inv = pinv(C_com);
C_mm = C_com_inv(1:k_com, 1:k_com);
C_nm = C_com_inv(1:k_com, end);

C_mm_inv = pinv(C_mm); 

scale_model.C_nm = C_nm;
scale_model.C_mm_inv = C_mm_inv;

% mu = -C_mm_inv*C_nm*[B_c; B_ring];
% [e_vec, e_val] = eig(C_mm_inv);
% B_tex = e_vec*sqrtm(e_val)*randn(k_tex, 1) + mu;