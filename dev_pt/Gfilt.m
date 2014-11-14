function [G,Gx,Gy,Gxx,Gxy,Gyy] = Gfilt(sigma,halfwidth)

% set debug flag
f_debug = false; if nargin==0, f_debug = true; end

if nargin<1, sigma = 9; end
if nargin<2, halfwidth = 3*sigma; end

sigmasq = sigma*sigma;
x	= (-halfwidth:halfwidth)';
g	= exp(-0.5* (x.*x)/sigmasq);
dg	= -x/sigmasq .* g;
ddg	= (-1/sigmasq * g) - (x/sigmasq .* dg);

G	= g*g';
Gx	= g*dg';
Gy	= dg*g';
Gxx	= g*ddg';
Gxy	= dg*dg';
Gyy	= ddg*g';

if ~f_debug, return; end

figure(1); clf; colormap(gray(256)); sbsz = [4,3];
subplot(3,3,1); imagesc(G); axis('image','xy');
% 	subplot(3,3,2); plot(x,g,'b',x,dg,'r',x,ddg,'g');
subplot(3,3,4); imagesc(Gx); axis('image','xy');
subplot(3,3,5); imagesc(Gy); axis('image','xy');
subplot(3,3,7); imagesc(Gxx); axis('image','xy');
subplot(3,3,8); imagesc(Gxy); axis('image','xy');
subplot(3,3,9); imagesc(Gyy); axis('image','xy');

t = rand*pi-pi/2
subplot(3,3,6); 
	imagesc(Gxx*cos(t)*cos(t)+Gyy*sin(t)*sin(t)+Gxy*sin(2*t)); axis('image','xy');
	c = halfwidth+1;
	hold on; plot(c+[0,3*sigma*cos(t)],c+[0,3*sigma*sin(t)],'r-');

[xx,yy] = meshgrid(x,x);
sigma = sigma/sqrt(2);
sigmasq = sigma*sigma;
x	= (-halfwidth:halfwidth)';
g2	= exp(-0.5* (x.*x)/sigmasq);
G2	= g2*g2';
figure(2); clf; colormap(gray(256)); sbsz = [7,2];
	mysubplot(sbsz); imagesc(xx.*yy); axis('image','xy');
	mysubplot(sbsz); imagesc(xx.*xx-yy.*yy); axis('image','xy');
	mysubplot(sbsz); imagesc(Gxy); axis('image','xy');
	mysubplot(sbsz); imagesc(Gxx-Gyy); axis('image','xy');
	mysubplot(sbsz); imagesc(G.*xx.*yy); axis('image','xy');
	mysubplot(sbsz); imagesc(G.*(xx.*xx-yy.*yy)); axis('image','xy');
	mysubplot(sbsz); imagesc(Gx.*Gy); axis('image','xy');
	mysubplot(sbsz); imagesc(Gx.*Gx-Gy.*Gy); axis('image','xy');
	mysubplot(sbsz); imagesc(Gx.*Gx); axis('image','xy');
	mysubplot(sbsz); imagesc(Gy.*Gy); axis('image','xy');
	tmp1 = conv2(Gx.*Gx,G2,'same');
	tmp2 = conv2(Gy.*Gy,G2,'same');
	mysubplot(sbsz); imagesc(tmp1); axis('image','xy');
	mysubplot(sbsz); imagesc(Gxx); axis('image','xy');
	mysubplot(sbsz); imagesc(tmp1-tmp2); axis('image','xy');
	mysubplot(sbsz); imagesc(Gxx-Gyy); axis('image','xy');
	
	

figpath = 'U:\projects\mammography\figs\filtering\';
imwrite(uint8(normim(Gx)*255),[figpath,'Gx.png']);
imwrite(uint8(normim(Gy)*255),[figpath,'Gy.png']);

imwrite(uint8(normim(Gxx)*255),[figpath,'Gxx.png']);
imwrite(uint8(normim(Gyy)*255),[figpath,'Gyy.png']);
imwrite(uint8(normim(Gxy)*255),[figpath,'Gxy.png']);
imwrite(uint8(normim(Gxx-Gyy)*255),[figpath,'Gxx-Gyy.png']);

imwrite(uint8(normim(Gx.*Gx)*255),[figpath,'GxGx.png']);
imwrite(uint8(normim(Gy.*Gy)*255),[figpath,'GyGy.png']);
imwrite(uint8(normim(Gx.*Gy)*255),[figpath,'GxGy.png']);
imwrite(uint8(normim(Gx.*Gx-Gy.*Gy)*255),[figpath,'GxGx-GyGy.png']);
	
clear G Gx Gy Gxx Gxy Gyy;

