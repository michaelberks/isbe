function flow_pyramid = build_flow_pyramid(flow, nLevels)

flow_pyramid = build_image_pyramid(flow, nLevels);

scl = 2;
for i = 2:nLevels
    flow_pyramid{i} = flow_pyramid{i} / scl;
    scl = scl * 2;
end

