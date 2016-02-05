function [xA,yA,xB,yB,xC,yC,xM,yM,c1n,c2n,lp]=brstpar(img)
%BRSTPAR pics the breast parameters manually from the image.
%   [XA,YA,XB,YB,XC,YC,XN,YN]=BRSTPAR(IMG) computes the breast parameters 
%   (XA,YA), (XB,YB), (XC,YC), (XN,YN) from the image IMG from for which 
%   reference points are given manually.  
%
%   [...,C1,C2,L]=BRSTPAR(IMG) additionally returns the breast boundary 
%   parabolae C1 and C2 and the pectoral line L. 

bAutoNormal = 0;

figure; clf;

imshow(img,[]);hold on;

disp('Please give the nipple point (A)');
[x,y]=ginput(1);
xA=x(1);
yA=y(1);
plot(xA,yA,'r+');

disp('Please give a point (U) on the upper breast boundary');
[x,y]=ginput(1);
xU=x(1);
yU=y(1);
plot(xU,yU,'+');

disp('Please give a point (L) on the lower breast boundary');
[x,y]=ginput(1);
xL=x(1);
yL=y(1);
plot(xL,yL,'+');

if ~bAutoNormal
    disp('Please give a breast point (M) on the normal direction from the nipple');
    [x,y]=ginput(1);
    xM=x(1);
    yM=y(1);
    plot(xM,yM,'r+');
end

disp('Please give three points on the pectoral muscle')
x=zeros(3,1); y=zeros(3,1);
for idx=1:3
  [xt,yt]=ginput(1);
  x(idx)=xt;
  y(idx)=yt;
  plot(x(idx),y(idx),'+')
end

pol=polyfit(x,y,1);
lp=[pol(1); -1; pol(2)]; %Pectoral line

% Find automatically the normal to the nipple point
if bAutoNormal
    %Intersection of the pectoral line perpendicular and the nipple line:
    xp=cross(cross([lp(1);lp(2);0],[xA;yA;1]),lp); 
    xp=reshape(xp/xp(3),3,1);
    
    %Breast boundary normal at the nipple
    n= xp(1:2) -[xA;yA]; 
    n=n/norm(n);    
    xM=xp(1);
    yM=xp(2);    
else
    n=[xM-xA;yM-yA];
    n=n/norm(n);    
end
%disp('Plotting the boundary parabolae...');
c1n=parabmin2([xA;xU],[yA;yU],n(1),n(2));
c2n=parabmin2([xA;xL],[yA;yL],n(1),n(2));

%plotconic(c1n);
%plotconic(c2n);

[x1,x2]=coniclinemeet(c1n,lp);
x1=x1/x1(3); x2=x2/x2(3);
laxis=cross([n;0],[xA;yA;1]);
if sign(laxis'*x1)==sign(laxis'*[xU;yU;1])
  xB=x1(1); yB=x1(2);
else
  xB=x2(1); yB=x2(2);  
end

[x1,x2]=coniclinemeet(c2n,lp);
x1=x1/x1(3); x2=x2/x2(3);
if sign(laxis'*x1)==sign(laxis'*[xL;yL;1])
  xC=x1(1); yC=x1(2);
else
  xC=x2(1); yC=x2(2);  
end

axis([0 ceil(max(max(size(img,2),xB),xC))...
      0 ceil(max(max(size(img,1),yB),yC))]);
plot(xB,yB,'r+');
plot(xC,yC,'r+');
plot([xB xC],[yB yC],'y-');

