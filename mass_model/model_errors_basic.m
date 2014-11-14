% MODEL_ERRORS_BASIC
%    [er_shape er_tex er_scale er_shape_scale B_shape_c B_tex_c B_scale_c]= 
%                   = model_errors_basic(mass_model, weights)
%
%    inputs:
%       mass_model  - 
%
%       optional:
%       weights     - 
%
%    outputs:
%
%    notes: % Formerly model_errors
%
% Copyright: (C) 2006-2008 University of Manchester
% Author: Michael Berks

function [er_shape er_tex er_scale er_shape_scale B_shape_c B_tex_c B_scale_c] =...
    model_errors_basic(mass_model, weights)

    if nargin < 2
        W_shape = mass_model.W_shape;
        W_tex   = mass_model.W_tex;
        W_scale = mass_model.W_scale;
        P_com = mass_model.P_com;
        B_com = mass_model.B_com;
    else
        W_shape = weights(1);
        W_tex   = weights(2);
        W_scale = 1;
        combined_data = [W_shape*mass_model.B_shape;...
            W_tex*mass_model.B_tex; W_scale*mass_model.B_scale]';
        [mean_com, P_com, B_com] = pca(combined_data, 0.98);
    end    
    
    k_shape     = size(mass_model.B_shape, 1);
    k_tex       = size(mass_model.B_tex, 1);
    
    Q_shape = P_com(1:k_shape,:); 
    Q_tex = P_com(k_shape+1:k_shape + k_tex,:);
    Q_scale = P_com(end, :);
    
    B_shape_c = inv(W_shape)*Q_shape*B_com;
    B_tex_c = inv(W_tex)*Q_tex*B_com;
    B_scale_c = Q_scale*B_com / W_scale;
    
    [N ssv] = size(mass_model.X_shape);
    ssv = ssv/2;
    
    X_shape_old = mass_model.P_shape*mass_model.B_shape + ...
        repmat(mass_model.mean_shape', 1, N);
    X_scale_old = mass_model.P_scale*mass_model.B_scale +...
        repmat(mass_model.mean_scale', 1, N);
    X_tex_old = mass_model.P_tex*mass_model.B_tex+...
        repmat(mass_model.mean_tex', 1, N);
    
    X_shape_new = mass_model.P_shape*B_shape_c + ...
        repmat(mass_model.mean_shape', 1, N);
    X_scale_new = mass_model.P_scale*B_scale_c +...
        repmat(mass_model.mean_scale', 1, N);
    X_tex_new = mass_model.P_tex*B_tex_c+...
        repmat(mass_model.mean_tex', 1, N);
    
    er_shape = X_shape_new - X_shape_old;
    er_shape = mean(sqrt(er_shape(1:ssv,:).^2 + er_shape(ssv+1:end,:).^2));
    er_tex = mean(abs(X_tex_new - X_tex_old));
    er_scale = X_scale_new - X_scale_old;
    er_shape_scale = repmat(X_scale_new, 2*ssv,1).*X_shape_new...
        - repmat(X_scale_old, 2*ssv, 1).*X_shape_old;
    er_shape_scale = mean(sqrt(er_shape_scale(1:ssv,:).^2 ...
        + er_shape_scale(ssv+1:end,:).^2));
    clear X_shape_new X_shape_old X_tex_new X_tex_old X_scale_new X_scale_old;
    

end