function [thickness_profile] = thickness_from_markers(x_pts, y_pts, M_mag, C_mag)
% Convert marker co-ordinates to breast thicknesses - we've changed this now
% so only relevant M_mag and C_mag values are passed in, no need for
% faffing
%--------------------------------------------------------------------------

distances = sqrt((x_pts(1,:)-x_pts(2,:)).^2 + (y_pts(1,:)-y_pts(2,:)).^2);
xpositions = mean(x_pts);

x_b = (distances - C_mag) ./M_mag;
thickness_profile = [xpositions' x_b'];
