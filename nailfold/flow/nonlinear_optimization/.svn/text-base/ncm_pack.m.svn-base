function [var_vec, sz_vec, lbound, ubound] = ncm_pack(variables)
% Pack structure of variables into a single vector that can be used within a
% nonlinear optimization routine.
% Variables that are not being used have a size of zero.

[vs, vt] = varnames();

nvs = length(vs);
nvt = length(vt);

sz_vec = zeros(nvs+nvt, 1);

% Add spatial variables such that those corresponding to the same pixel are
% consecutive. This gives a better block diagonal structure to the 
% (approximated) Hessian.
vs_mat = [];
lb_mat = [];
ub_mat = [];

row = 0;
for i = 1:nvs
    if isfield(variables, vs{i})
        row = row + 1;

        var_value = variables.(vs{i});
        n = numel(var_value);

        [lb, ub] = varlimits(vs{i});
        lb_mat(row, :) = lb*ones(1, numel(var_value));
        ub_mat(row, :) = ub*ones(1, numel(var_value));

        if varsquash(vs{i})
            % Convert a limited variable (e.g. occupancy) to an unlimited
            % one
            var_value(:) = var_value(:) - lb;
            var_value(:) = var_value(:) / (ub-lb);
            vs_mat(row, :) = invsigmoid(var_value(:)');
        else
            vs_mat(row, :) = var_value(:)';
        end
        
        sz_vec(i) = n; 
    end
end

var_vec = vs_mat(:);
lbound = lb_mat(:);
ubound = ub_mat(:);

% Add temporal variables on last.
for i = 1:nvt
    if isfield(variables, vt{i})
        var_value = variables.(vt{i});
        n = numel(var_value);
        
        [lb, ub] = varlimits(vt{i});
        lbound(end+1:end+n) = lb;
        ubound(end+1:end+n) = ub;

        if varsquash(vt{i})
            % Convert a limited variable (e.g. occupancy) to an unlimited
            % one
            var_value(:) = var_value(:) - lb;
            var_value(:) = var_value(:) / (ub-lb);
            var_vec(end+1:end+n) = invsigmoid(var_value(:));
        else
            var_vec(end+1:end+n) = var_value(:);
        end

        sz_vec(nvs+i) = n;
    end
end
