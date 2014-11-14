function c = vessel_curvature(vessel,step)

if nargin<2, step = 1; end

theta = nan(size(vessel,1),1);
for i = 1+step:size(vessel,1)-step
	v1 = vessel(i,:)-vessel(i-step,:);
	v2 = vessel(i+step,:)-vessel(i,:);
	theta(i) = acos((v1*v2') / (norm(v1)*norm(v2)));
	if isnan(theta(i)), keyboard; end
end

c = [diff(theta); nan];

