function [max_icp] = compute_max_icp(icp_level)

[max_icp] = max(icp_level, [], 3);

icp_abs = zeros(size(max_icp));
icp_phase = angle(max_icp);
for pp = 1:6
    icp_abs = ...
        icp_abs + ...
        abs(abs(icp_level(:,:,pp)) .* ...
            cos(angle(icp_level(:,:,pp)) - icp_phase));
end
max_icp = icp_abs .* exp(i*icp_phase);

    