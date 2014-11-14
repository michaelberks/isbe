function test_flow_scale(variables, observations)
% Exhaustive testing of flow rate

f0 = variables.f;
evec = [];
scales = linspace(-1,3,21);
for s = scales
    variables.f = s*f0;
    [v, sz_vec] = pack(variables);
    evec(end+1) = error_val(v, sz_vec, observations);
end
figure(3); clf;
    plot(scales, evec);