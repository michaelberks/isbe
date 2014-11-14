function flowmap = flowmap_from(p, base_flow)

if isfield(p, 'uu')
    flowmap = complex(p.uu, ...
                      p.vv);
else
    flowmap = complex(p.F .* cos(p.theta), ...
                      p.F .* sin(p.theta));
end

if ~isempty(base_flow) && ...
   all(size(flowmap) == size(base_flow))
    flowmap = base_flow + flowmap;
end
