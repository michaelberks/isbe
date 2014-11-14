% COMBINE_MODEL - recompute combined appearance model using new set of weights 
%    [mass_model] = combine_model(mass_model, weights, file_out)
%
%    inputs:
%       mass_model  - structure containing appearance model
%       weights     - 3x1 vector of new shape, texture and scale weights
%       file_out    - file path to save new model to
%
%    outputs:
%       mass_model - structure containing recomputed appearance model
%
%    notes:
%
% Copyright: (C) 2006-2008 University of Manchester
% Author: Michael Berks

function mass_model = combine_model(mass_model, weights, file_out)

    switch length(weights)
        case 2
            W_shape = weights(1);
            W_tex   = weights(2);
            W_scale = 1;
        case 3
            W_shape = weights(1);
            W_tex   = weights(2);
            W_scale = weights(3);
        case size(mass_model.combined_data, 2)
            W_shape = diag(weights(1:length(mass_model.L_shape)));
            W_tex = diag(weights(length(mass_model.L_shape)+1:end-1));
            W_scale = weights(end);
        otherwise
            error('Weight vector is of incorrect length');
    end

    mass_model.W_shape = W_shape;
    mass_model.W_tex = W_tex;
    mass_model.W_scale = W_scale;
    
    combined_data = [W_shape*mass_model.B_shape; W_tex*mass_model.B_tex; W_scale*mass_model.B_scale]';
    mass_model.combined_data = combined_data;

    [mean_com, P_com, B_com, L_com] = pca(combined_data, 0.98);

    mass_model.mean_com = mean_com;
    mass_model.P_com = P_com;
    mass_model.B_com = B_com;
    mass_model.L_com = L_com;

    mass_model.progress = '7: Combined model complete, function successful!';
    
    if nargin > 2
        save(file_out, 'mass_model');
    end

end