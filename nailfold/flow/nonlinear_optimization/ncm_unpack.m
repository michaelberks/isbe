function p = ncm_unpack(var_vec, sz_vec)
% Split a packed vector into a structure of fields.
% Variables with a size of zero are not extracted.

p = [];

[vs, vt] = varnames();

nvs = length(vs);
nvt = length(vt);

vs_sz = sz_vec(1:nvs);
vs_sz = vs_sz(vs_sz ~= 0);

if ~isempty(vs_sz)
    if any(vs_sz ~= vs_sz(1))
        error('Spatial variables must have same size');
    end
    
    nvs_nz = numel(vs_sz);
    vs_mat = reshape(var_vec(1:sum(vs_sz)), [nvs_nz, vs_sz(1)]);

    row = 1;
    for i = 1:nvs
        if (sz_vec(i) > 0)
            if varsquash(vs{i})
                [lb, ub] = varlimits(vs{i});            
                p.(vs{i}) = lb + (ub-lb)*sigmoid(vs_mat(row, :)');
            else
                p.(vs{i}) = vs_mat(row, :)';
            end
            row = row + 1;
        end
    end
end

first = sum(vs_sz)+1;
for i = 1:nvt
    if (sz_vec(nvs+i) > 0)
        last = first + sz_vec(nvs+i) - 1;
        if varsquash(vt{i})
            [lb, ub] = varlimits(vt{i});
            p.(vt{i}) = lb + (ub-lb)*sigmoid(var_vec(first:last));
        else
            p.(vt{i}) = var_vec(first:last);
        end
        first = last + 1;
    end
end
