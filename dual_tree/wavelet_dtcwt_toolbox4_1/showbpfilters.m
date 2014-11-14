% showbpfilters.m
%
%Show the 2-D responses of the unmodified and modified DT CWT

nlevels = 4;
m = 5;
xs = [2,6]*m*(2^(nlevels)); % Size of image
x = zeros(xs);
% x = rand(xs);

% x(xs/2) = 1;

[Yl,Yh] = dtwavexfm2(x,nlevels,biort,qshift);

for k = 1:6,
    a = (m+1)/2;
    b = a + (k-1)*m;
    Yh{nlevels}(a,b,k) = 100;
    Yh{nlevels}(a+m,b,k) = 100*j;
end

Z = dtwaveifm2(Yl,Yh,biort,qshift); 
setfig(11);
image(uint8(Z*20+128)); colormap(gray(256))
axis image
axis off
text(0.5*xs(2),-0.1*xs(1),'(a) Dual-Tree Complex Wavelets: Real Part','horiz','c');
text(0.5*xs(2),1.1*xs(1),'Imaginary Part','horiz','c');
label = [' 15 ';' 45 ';' 75 ';'105 ';'135 ';'165 '];
for k = 1:6
    text(xs(2)*(k-0.5)/6,0.5*xs(1),label(k,:),'horiz','c');
end

print -depsc impresp2d.eps

% sbcor = j*[1 -0.9115 1 j -0.9115*j j]; % Phase and amplitude correction for individual subbands.
sbcor = [j -0.9115*j j -1 0.9115 -1]; % Phase and amplitude correction for individual subbands.
for k = 1:6, Yh{nlevels}(:,:,k) = Yh{nlevels}(:,:,k) / sbcor(k); end

Zb = dtwaveifm2b(Yl,Yh,biort,[qshift '_bp']); 
setfig(12);
image(uint8(Zb*20+128)); colormap(gray(256))
axis image
axis off
text(0.5*xs(2),-0.1*xs(1),'(b) Modified Complex Wavelets: Real Part','horiz','c');
text(0.5*xs(2),1.1*xs(1),'Imaginary Part','horiz','c');
label = [' 15 ';' 45 ';' 75 ';'105 ';'135 ';'165 '];
for k = 1:6
    text(xs(2)*(k-0.5)/6,0.5*xs(1),label(k,:),'horiz','c');
end

print -depsc impresp2dbp.eps

setfig(13)
f = [-128:127]*4/256;
plot(f,abs(fftshift(fft(hc,256))),'--r','Linewidth',2); 
hold on
plot(f,abs(fftshift(fft(hcnew,256))),'-b','Linewidth',2); 
hold off
grid on
set(gca,'pos',[0.1550 0.300 0.7250 0.35],'Fontsize',14)
title('(c) Frequency responses of original and modified 1-D filters');
xlabel('frequency / output sample rate');
text(1.1,1.98,'original','horiz','l','vert','t');
text(0.25,1.98,'modified','horiz','r','vert','t');

print -depsc freqresp1d.eps


