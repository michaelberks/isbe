xmax = 12;
x = linspace(-xmax,xmax,601);
dx = x(2)-x(1);

y0 = zeros(size(x));
y0(abs(x)<1) = 1;
y0 = y0/(sum(y0)*dx);

g = exp(-0.5*x.*x/1); g = g/(sum(g)*dx);

y = y0;

figure(1); clf; hold on;
	plot(x,y,'b-');
	plot(x,g,'r-');
	axis([-xmax,xmax,0,0.6]);

for i = 2:30
	i
	y = conv2(y,y0,'same')*dx;
	
	s = 0.55*log(i-1)/log(2);
	g = exp(-0.5*x.*x/(s*s));	g = g/(sum(g)*dx);
	
	figure(1); clf; hold on;
	plot(x,y,'b-');
	plot(x,g,'r-');
	axis([-xmax,xmax,0,0.6]);
	pause(0.04);
end