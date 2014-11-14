clc;

figpath = 'U:\projects\mammography\figs\filtering\';

%% 2D
sigma = 2;
[ddG,G,dG] = Wfilt(0,sigma);
fhw = (size(G,1)-1)/2; % filter halfwidth

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
[xx,yy] = meshgrid(-N:N,-N:N);
r		= sqrt(xx.*xx+yy.*yy);
T		= atan2(yy,xx);
imgtype = 'fingerprint';
switch imgtype
	case 'square',
		img = -ones(size(xx));
		img(abs(xx)<=(N*0.75) & abs(yy)<=(N*0.75)) = 1;
		img(abs(xx)<=(N*0.7) & abs(yy)<=(N*0.7)) = -1;
		origt = pi/4+sign(xx.*xx-yy.*yy)*pi/4;
	case 'disk',
		img = -ones(size(xx));
		img(r<=(N*0.75)) = 1;
		img(r<=(N*0.7)) = -1;
		origt = mod(atan2(yy,xx) + pi/2,pi);
	case 'horizline',
		img = -ones(size(xx));
		img(abs(yy)<3) = 1;
		origt = zeros(size(img));
	case 'vertline',
		img = -ones(size(xx));
		img(abs(xx)<3) = 1;
		origt = pi/2*ones(size(img));
	case 'karstest',
		origt = kars_test_image;
		img = -ones(size(origt)); img(origt>0) = 1;
		N = round(size(img,1)/2);
	case 'fingerprint',
		img = double(imread('U:\data\fingerprints\neurotechnology\012_3_1.png'));
		img = double(sign(192-img));
		origt = zeros(size(img));
end

% centre pixel
c = round(numel(img)/2);

% create mask of pixels to keep
mask = (img>0);

% invert?
% img = -img;

img0 = img;
% img = conv2(img,G,'same');
img = normim(img)*2 - 1;
img = img + 0.0*randn(size(img));

dummyval = nan;
origt(~mask) = dummyval;
% origt(mask) = pi;

	mysubplot(sbsz); imagesc(img); axis('image','xy');
	dimg = conv2(img,dG1,'same');
	mysubplot(sbsz); imagesc(dimg); axis('image','xy');
	ddimg = conv2(img,W1,'same');
	mysubplot(sbsz); imagesc(ddimg); axis('image','xy');
	
%% equivalence between three ddG filters and others

% get responses to various filters
[G,Gx,Gy,Gxx,Gxy,Gyy] = Gfilt(sigma);
Ix = conv2(img,Gx,'same');
Iy = conv2(img,Gy,'same');
Ixx = conv2(img,Gxx,'same');
Ixy = conv2(img,Gxy,'same');
Iyy = conv2(img,Gyy,'same');
Iw0 = conv2(img,W0,'same');
Iw1 = conv2(img,W1,'same');
Iw2 = conv2(img,W2,'same');

% use first derivatives without squared gradients
KcosT = Ix; KsinT = Iy;
ori3 = mod(atan2(KsinT,KcosT) + pi/2,pi);
ori3(~mask) = dummyval;
err3 = mean(abs(ori_error(origt(mask),ori3(mask))));

% use first derivatives with squared gradients
Ksin2Tt2 = 2*Ix.*Iy; Kcos2Tt2 = Ix.*Ix - Iy.*Iy;
ori3b = 0.5*(atan2(Ksin2Tt2,Kcos2Tt2)) + pi/2;
ori3b(~mask) = dummyval;
err3b = mean(abs(ori_error(origt(mask),ori3b(mask))));

% apply Kars to circle image
Ksin2T = sqrt(3)*(Iw2-Iw1); Kcos2T = Iw1+Iw2-2*Iw0;
ori1 = 0.5*(atan2(Ksin2T,Kcos2T)) + pi/2;
% Wori1 = (1-cos(2*ori1))*Ksin2T - sin(2*ori1)*Kcos2T + 3*FI0;
% Wori1perp = (1-cos(2*ori1+pi))*Ksin2T - sin(2*ori1+pi)*Kcos2T + 3*FI0;
% ori1(Wori1<Wori1perp) = ori1(Wori1<Wori1perp) + pi/2;
% ori1 = mod(ori1,pi); % for display only

ori1(~mask) = dummyval;
err1 = mean(abs(ori_error(origt(mask),ori1(mask))));

% try same thing with equivalent cloverleaf filters
Ksin2Tb = 2*Ixy; Kcos2Tb = (Ixx-Iyy);
ori2 = 0.5*(atan(Ksin2Tb./Kcos2Tb));
R	= Ixx.*cos(ori2).*cos(ori2) + Iyy.*sin(ori2).*sin(ori2) + Ixy.*sin(2*ori2);
Rp	= Ixx + Iyy - R; % response at perpendicular
imin = find(abs(Rp)>abs(R)); % minima of negative responses
ori2(imin) = ori2(imin) + pi/2;
ori2 = ori2 + pi/2; % orientation is perpendicular to gradient
ori2(abs(R-Rp)<1e-6) = dummyval; % exclude regions that are completely flat
ori2(~mask) = dummyval;
err2 = mean(abs(ori_error(origt(mask),ori2(mask))));
ori2 = mod(ori2,pi); % for display only

% try again with Haar-like approximation
[fx,fy] = meshgrid(-fhw:fhw,-fhw:fhw);
h_mask	= abs(fx)<=fhw/sqrt(2) & abs(fy)<=fhw/sqrt(2);
f1e		= sign(fx.*fy) .* h_mask;
h_mask	= abs(fx+fy)<=fhw & abs(fx-fy)<=fhw;
f2e		= sign(fx.*fx-fy.*fy) .* h_mask;
Ksin2Th = conv2(img,f1e,'same');
Kcos2Th = conv2(img,f2e,'same');
ori2h = 0.5*( atan2(Ksin2Th,Kcos2Th) ) ;
ori2h(~mask) = dummyval;
err2h = mean(abs(ori_error(origt(mask),ori2h(mask))));
ori2h = mod(ori2h,pi); % for display only

[	err3 err3b ... % 1st derivs
	err1 err2  ... % 2nd derivs
	err2h ... % haar approx
]

[	KcosT(c) KsinT(c);
	Kcos2Tt2(c) Ksin2Tt2(c);
	Kcos2T(c) Ksin2T(c);
	Kcos2Tb(c) Ksin2Tb(c);
	Kcos2Th(c) Ksin2Th(c);]

figure(2); clf; sbsz = [2,2]; colormap(gray(256));
	mysubplot(sbsz); imagesc(Gxy); axis('image','xy');
	mysubplot(sbsz); imagesc(Gxx-Gyy); axis('image','xy');
	mysubplot(sbsz); imagesc(f1e); axis('image','xy');
	mysubplot(sbsz); imagesc(f2e); axis('image','xy');

cmap = hsv(255); cmap = [0 0 0; cmap];
figure(3); clf; sbsz = [2,3]; colormap(cmap);
	% original
	mysubplot(sbsz); image(uint8(1+255*origt/pi)); axis('image','xy'); title('GT');
	% first derivatives (without & with squaring)
	mysubplot(sbsz); image(uint8(1+255*ori3/pi)); axis('image','xy'); title('1st unsquared');
	mysubplot(sbsz); image(uint8(1+255*ori3b/pi)); axis('image','xy'); title('1st squared');
	% karssemeijer's method
	mysubplot(sbsz); image(uint8(1+255*ori1/pi)); axis('image','xy'); title('kars');
	% equivalent 'cloverleaf' filters
	mysubplot(sbsz); image(uint8(1+255*ori2/pi)); axis('image','xy'); title('cloverleaf');
	% Haar-like approximations
	mysubplot(sbsz); image(uint8(1+255*ori2h/pi)); axis('image','xy'); title('Haar');

	
% output images
imwrite(uint8(normim(W0)*255),[figpath,'W0.png']);
imwrite(uint8(normim(W1)*255),[figpath,'W1.png']);
imwrite(uint8(normim(W2)*255),[figpath,'W2.png']);
imwrite(uint8(normim(Gx)*255),[figpath,'Gx.png']);
imwrite(uint8(normim(Gy)*255),[figpath,'Gy.png']);
imwrite(uint8(normim(Gxx)*255),[figpath,'Gxx.png']);
imwrite(uint8(normim(Gyy)*255),[figpath,'Gyy.png']);
imwrite(uint8(normim(Gxy)*255),[figpath,'Gxy.png']);
imwrite(uint8(normim(Gxx-Gyy)*255),[figpath,'Gxx-Gyy.png']);
imwrite(uint8(normim(f1e)*255),[figpath,'haar_sin.png']);
imwrite(uint8(normim(f2e)*255),[figpath,'haar_cos.png']);

% outputs from orientation estimation
if ~exist([figpath,imgtype],'dir')
	mkdir([figpath,imgtype]);
end
imwrite(uint8(normim(img)*255),[figpath,imgtype,'\input.png']);
imwrite(uint8(1+255*origt/pi),cmap,[figpath,imgtype,'\ori_gt.png']);
imwrite(uint8(1+255*ori3/pi),cmap,[figpath,imgtype,'\ori_1st.png']);
imwrite(uint8(1+255*ori3b/pi),cmap,[figpath,imgtype,'\ori_1st_sqr.png']);
imwrite(uint8(1+255*ori1/pi),cmap,[figpath,imgtype,'\ori_kars.png']);
imwrite(uint8(1+255*ori2/pi),cmap,[figpath,imgtype,'\ori_clover.png']);
imwrite(uint8(1+255*ori2h/pi),cmap,[figpath,imgtype,'\ori_haar.png']);

	