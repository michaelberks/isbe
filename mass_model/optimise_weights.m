% OPTIMISE_WEIGHTS
%    [weights_out, output] =... 
%       optimise_weights(mass_files, weights_in, model_path, mass_path)
%
%    inputs:
%       mass_files  - 
%       weights_in  -
%       model_path  -
%       mass_path   -
%
%    outputs:
%       weights_out -
%       output      -
%    notes: Look at using MDL as objective function?
%
% Copyright: (C) 2006-2008 University of Manchester
% Author: Michael Berks

function [weights_out, output] =...
    optimise_weights(mass_model, mass_files, weights_in, save_file, indices)    
      
    options = optimset('LargeScale', 'off',...
        'MaxIter', 1000, 'MaxFunEvals', 5000);
    [weights_out, fval, exitflag, output] =...
        fmincon(@myfun,weights_in,[],[],[],[], [1e-6 1e-6], [inf inf],...
            [], options);
    %weights_out = fminunc(@myfun,weights_in);
    
    function er = myfun(x)
%         [e_c] = model_errors_loo(mass_files, x, model_path, indices, mass_path);
        [e_c] = ...
            model_errors_loo2(mass_model, mass_files, 'weights', x, 'indices', indices);
        er = mean(e_c.weights);
        save(save_file, 'x');
    end
    save(save_file, 'weights_out', 'output');
end