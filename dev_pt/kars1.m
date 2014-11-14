figpath = 'U:\projects\mammography\figs\filtering\';

%% 1D
% x = linspace(-3,3,101);
% s = 1;
% gx = (2*pi*s*s)^-0.5 * exp(-x.*x/(2*s*s));
% dgx = -1/(s*s)*x.*gx;
% ddgx = -1/(s*s)*(-x.*x/(s*s) + 1) .* gx;
% figure(2); plot(x,gx,'r-',x,dgx,'g-',x,ddgx,'b-');
% return

%% 2D
sigma = 7;
[ddG,G,dG] = Wfilt(rand*pi,sigma);
[ddG,G,dG] = Wfilt(0,sigma);

sbsz = [3,3];
figure(1); clf; colormap(gray(256));
	mysubplot(sbsz); imagesc(G); axis('equal','tight');
	mysubplot(sbsz); imagesc(dG); axis('equal','tight'); 
	mysubplot(sbsz); imagesc(ddG); axis('equal','tight');
	
[W0,G0,dG0] = Wfilt(0,sigma);
[W1,G1,dG1] = Wfilt(pi/3,sigma);
[W2,G2,dG2] = Wfilt(2*pi/3,sigma);
	
	mysubplot(sbsz); imagesc(W0); axis('equal','tight');
	mysubplot(sbsz); imagesc(W1); axis('equal','tight');
	mysubplot(sbsz); imagesc(W2); axis('equal','tight');
	
% generate an image of a disk
N = 100;
[xx,yy] = meshgrid(-N:N,N:-1:-N);
img = ones(size(xx));
% r = sqrt(xx.*xx+yy.*yy);
% img(r<=(N*0.75)) = 1;
% img(r<=(N*0.7)) = -1;
% origt = mod(atan2(yy,xx) + pi/2,pi);

% img(abs(yy)<3) = 1;
% origt = zeros(size(img));

img(abs(xx)<3) = -1;
origt = pi/2*ones(size(img));

% origt = kars_test_image;
% img = -ones(size(origt)); img(origt>0) = 1;

img0 = img;
img = conv2(img,G,'same');
img = img-min(img(:));
img = img/max(img(:))*2 - 1;
img = img + 0.0*randn(size(img));

origt(img0<0) = NaN;

	mysubplot(sbsz); imagesc(img); axis('image');
	dimg = conv2(img,dG1,'same');
	mysubplot(sbsz); imagesc(dimg); axis('image');
	ddimg = conv2(img,W1,'same');
	mysubplot(sbsz); imagesc(ddimg); axis('image');
	
%% equivalence between three ddG filters and... a four-leafed clover???
f1	= sqrt(3)*(W2 - W1);
f2	= W1 + W2 - 2*W0;
	
[rank(W0) rank(W1) rank(W2)]
[rank(f1) rank(f2)]

c = sub2ind(size(img),N-1,N-1);

% apply Kars to circle image
FI0 = conv2(img,W0,'same');
FI1 = conv2(img,W1,'same');
FI2 = conv2(img,W2,'same');
Ksin2T = sqrt(3)*(FI2-FI1);
Kcos2T = FI1+FI2-2*FI0;
ori1 = 0.5*( atan2(Ksin2T,Kcos2T) ) + pi/2;
ori1(img0<0) = NaN;
err1 = mean(abs(ori_error(origt(img0>0),ori1(img0>0))));

% try same thing with equivalent filters
Ksin2Tb = conv2(img,f1,'same');
Kcos2Tb = conv2(img,f2,'same');
ori2 = 0.5*( atan2(Ksin2Tb,Kcos2Tb) ) + pi/2;
ori2(img0<0) = NaN;
err2 = mean(abs(ori_error(origt(img0>0),ori2(img0>0))));

% try again with Haar-like approximation
f1e	= sign(f1); f1e(:,(size(f1,1)+1)/2) = 0; f1e((size(f1,2)+1)/2,:) = 0;
f2e = sign(f2); f2e(1:size(f2e,1)+1:end) = 0; f2e(size(f2e,1):size(f2e,1)-1:end) = 0;
Ksin2Th = conv2(img,f1e,'same');
Kcos2Th = conv2(img,f2e,'same');
ori2h = 0.5*( atan2(Ksin2Th,Kcos2Th) ) + pi/2;
ori2h(img0<0) = NaN;
err2h = mean(abs(ori_error(origt(img0>0),ori2h(img0>0))));

% try again with smoothed product approximation
f1d = dG(9,:)'*dG(9,:); scl = f1./f1d;
scl = mean(scl(~isnan(scl) & ~isinf(abs(scl))));
f1d = f1d*scl;
Ksin2Td = conv2(img,f1d,'same');
Kcos2Td = conv2(img,f2,'same');
ori2d = 0.5*( atan2(Ksin2Td,Kcos2Td) ) + pi/2;
ori2d(img0<0) = NaN;
err2d = mean(abs(ori_error(origt(img0>0),ori2d(img0>0))));

% do individually
RcosT = conv2(img,dG,'same');
RsinT = conv2(img,dG','same');
Ksin2Tt = 2*RsinT.*RcosT;
Kcos2Tt = RcosT.*RcosT - RsinT.*RsinT;
ori3 = 0.5*( atan2(Ksin2Tt,Kcos2Tt) ) + pi/2;
ori3(img0<0) = NaN;
err3 = mean(abs(ori_error(origt(img0>0),ori3(img0>0))));

[err1 err2 err2h err2d err3]

[	Ksin2T(c) Kcos2T(c);
	Ksin2Tb(c) Kcos2Tb(c);
	Ksin2Th(c) Kcos2Th(c);
	Ksin2Td(c) Kcos2Td(c);
	Ksin2Tt(c) Kcos2Tt(c) ]


figure(2); clf; sbsz = [2,2]; colormap(gray(256));
	mysubplot(sbsz); imagesc(f1); axis('image');
	mysubplot(sbsz); imagesc(f2); axis('image');
	mysubplot(sbsz); imagesc(f1d); axis('image');
	mysubplot(sbsz); imagesc(f2e); axis('image');

cmap = [0 0 0; hsv(255)];
figure(3); clf; sbsz = [2,3]; colormap(cmap);
	mysubplot(sbsz); imagesc(origt); axis('image');
	mysubplot(sbsz); imagesc(ori1); axis('image');
	mysubplot(sbsz); imagesc(ori2); axis('image');
	mysubplot(sbsz); imagesc(ori2h); axis('image');
	mysubplot(sbsz); imagesc(ori2d); axis('image');
	mysubplot(sbsz); imagesc(ori3); axis('image');

return
	
% output images
imwrite(uint8(normim(W0)*255),[figpath,'W0.png']);
imwrite(uint8(normim(W1)*255),[figpath,'W1.png']);
imwrite(uint8(normim(W2)*255),[figpath,'W2.png']);
imwrite(uint8(normim(f1)*255),[figpath,'f1.png']);
imwrite(uint8(normim(f2)*255),[figpath,'f2.png']);

% decompose f1 and f2 into separable filters
[U,S,V] = svd(f1,'econ');
for i = 1:rank(f1)
	f1sep = U(:,i)*S(i,i)*V(:,i)';
	filename = [figpath,sprintf('f1_%i.png',i)];
	imwrite(uint8(normim(f1sep)*255),filename);
end

[U,S,V] = svd(f2,'econ');
for i = 1:rank(f2)
	f2sep = U(:,i)*S(i,i)*V(:,i)';
	filename = [figpath,sprintf('f2_%i.png',i)];
	imwrite(uint8(normim(f2sep)*255),filename);
end

	