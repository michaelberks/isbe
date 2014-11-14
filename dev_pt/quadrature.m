fsz = 31;
g	= exp(-0.5*linspace(-3,3,fsz).^2)';
x	= 4*linspace(-pi,pi,31);
f1	= g*(cos(x).*g');
f2	= g*(sin(x).*g');

figure(2); clf; colormap(gray(256));
	subplot(2,4,2); plot(x,cos(x).*g'); axis('square',[-4*pi,4*pi,-1,1]);
	subplot(2,4,5); plot(1-g,x); axis('square',[0,1,-4*pi,4*pi]);
	subplot(2,4,6); image(uint8(255 * (f1+1)/2)); axis('image');
	subplot(2,4,3); plot(x,sin(x).*g'); axis('square',[-4*pi,4*pi,-1,1]);
	subplot(2,4,8); plot(g,x); axis('square',[0,1,-4*pi,4*pi]);
	subplot(2,4,7); image(uint8(255 * (f2+1)/2)); axis('image');
	
% 2D description as in Freeman and Adelson
W	= g*g';
Px1	= ones(fsz,1)*cos(x); f1 = W.*Px1;
Px2	= ones(fsz,1)*sin(x); f2 = W.*Px2;

for t = 1:360
	f = cosd(t)*f1 + sind(t)*f2;
	figure(1); clf; colormap(gray(256));
		image(uint8(255 * (f+1)/2)); axis('image');
		drawnow;
end
