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

    W_shape = weights(1);
    W_tex   = weights(2);
    W_scale = weights(3);

    mass_model.W_shape = W_shape;
    mass_model.W_tex = W_tex;
    mass_model.W_scale = W_scale;
    
    combined_data = [W_shape*mass_model.B_shape; W_tex*mass_model.B_tex; W_scale*mass_model.B_scale]';
    mass_model.combined_data = combined_data;

    [mean_c, P_c, B_c, L_c] = pca(combined_data, 0.98);

    mass_model.mean_c = mean_c;
    mass_model.P_c = P_c;
    mass_model.B_c = B_c;
    mass_model.L_c = L_c;

    mass_model.progress = '7: Combined model complete, function successful!';
    
    if nargin > 2
        save(file_out, 'mass_model');
    end

end