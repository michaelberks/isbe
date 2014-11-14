T = linspace(0,360,501); T(end) = [];
v = zeros(size(T));

for iT = 1:length(T)
	t = linspace(0,T(iT),100);
	v(iT) = abs(mean(complex(cosd(t),sind(t))));
end
v = 1-acos(v)/(pi/2);

[ignore,i(1)] = min(abs(v-0.5));
[ignore,i(2)] = min(abs(v-0.9));
[ignore,i(3)] = min(abs(v-0.99));
[ignore,i(4)] = min(abs(v-0.999));

T(i)

figure(1); clf; hold on;
plot(T,v,'b-');
plot([T(i); T(i)],1e6*[-1,1],'-','color',0.8*[1,1,1]);
ylim([0,1]);

graph(1); exportfig([asymmetryroot,'figs/dispersions']);