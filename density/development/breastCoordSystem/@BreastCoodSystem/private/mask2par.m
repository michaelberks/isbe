function [coordSystem]=mask2par(objCS, mask,xAinit,yAinit,lp)
%MASK2PAR computes the breast parameters from the mask
%   [XA,YA,XB,YB,XC,YC,C1,C2]=MASK2PAR(MASK,XA_INIT,YA_INIT,L) computes the breast
%   boundary approximation form MASK and the nipple approximation
%   (XA_INIT,YA_INIT). The nipple point will be the vertex of the
%   parabolae which will locate on the breast normal line, perpendicular to the
%   pectoral line.  
%
%   [...,CPOINTS]=MASK2CONIC(MASK,XA_INIT,YA_INIT,L) additionally
%   returns points on the boundary approximation.

lp = reshape(lp,length(lp),1);

%Intersection of the pectoral line perpendicular and the nipple line:
xp=cross(cross([lp(1);lp(2);0],[xAinit;yAinit;1]),lp); 
xp=reshape(xp/xp(3),3,1);

%First extract the breast boundary points
%Assuming that the boundary points of the labels '0' and '1' are taken
[ y, x] = extractCroppedBoundary(mask,objCS.labelBreast, objCS.labelBackground); % ,[xp(2) xp(1), [yAinit xAinit]

%Breast boundary normal at the nipple
n= xp(1:2) -[xAinit;yAinit]; 
n=n/norm(n);

R=[[n(2);-n(1)] n]'; %Rotation that puts the normal towards the y-axis

T=[R -R*[xAinit;yAinit]; 0 0 1]; %Transformation

Xt=T*[x';y';ones(size(x'))]; %Transformed boundary points
%lt=inv(T)'*l; %Transformed pectoral line


Xpt=T*xp; %Pectoral line point at the new y-axis
D=Xpt(2); %The distance of the pectoral line from the nipple

%The least squares fitting of the boundary points
%The parabolae are y=a_left*x+c and y=a_right*x+c 

%Define the weighting of the points
W=(D-Xt(2,:))/D;

 W(W<0) = 0;
 W(W>1) = 0;


idxL=find(Xt(1,:)<=0); %Left branch
idxR=find(Xt(1,:)>=0); %Right branch

%Measurement matrix
M=zeros(length(idxL)+length(idxR),3);
M(:,3)=[W(idxL).^(1/2)';W(idxR).^(1/2)']; %Homogenization
M(1:length(idxL),1)=W(idxL).^(1/2).*Xt(1,idxL).^2;
M(length(idxL)+1:end,2)=W(idxR).^(1/2).*Xt(1,idxR).^2;
Y=[W(idxL).^(1/2)'.*Xt(2,idxL)';W(idxR).^(1/2)'.*Xt(2,idxR)'];

A=pinv(M)*Y;
a_left=A(1); a_right=A(2); c=A(3);

%Transformation of the conics to the original coordinate system
CtL=[a_left 0 0; 0 0 -0.5; 0 -0.5 c];
CtR=[a_right 0 0; 0 0 -0.5; 0 -0.5 c];

if n(1)<0 %Check the orientation of the breast
  C1=T'*CtL*T;
  C2=T'*CtR*T;
else
  C1=T'*CtR*T;
  C2=T'*CtL*T;
end

%Transformation of the refined nipple position
XA=T\[0;c;1];
% XA=T\[0;1;1];  % JM
xA=XA(1); yA=XA(2);

%Finding the intersection points B and C on the pectoral line
[x1,x2]=coniclinemeet(C1,lp);
x1=x1/x1(3); x2=x2/x2(3);
if xp(2) >= x1(2) %Compare the y-coordinates to find the right solution
  xB=x1(1); yB=x1(2);
else
  xB=x2(1); yB=x2(2);  
end

[x1,x2]=coniclinemeet(C2,lp);
x1=x1/x1(3); x2=x2/x2(3);
if xp(2) <= x1(2)
  xC=x1(1); yC=x1(2);
else
  xC=x2(1); yC=x2(2);  
end

%Transformation of the points
[tmp,minIdx]=min(Xt(1,:));
[tmp,maxIdx]=max(Xt(1,:));

xleft=Xt(1,minIdx):0;
xright=0:Xt(1,maxIdx);
yleft=a_left*xleft.^2+c;
yright=a_right*xright.^2+c;
Cpoints=T\[xleft xright;yleft yright;ones(1,length(xleft)+length(xright))];

coordSystem = struct;
coordSystem.xA = xA;
coordSystem.yA = yA;
coordSystem.xB = xB;
coordSystem.yB = yB;           
coordSystem.xC = xC;
coordSystem.yC = yC;           
coordSystem.c1n = C1;
coordSystem.c2n = C2;
coordSystem.Cpoints = Cpoints;
end


