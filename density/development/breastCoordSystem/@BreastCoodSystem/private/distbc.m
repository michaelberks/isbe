function result=distbc(phi,X,Y,xA,yA,xB,yB,xC,yC,phi_0,d,c1,c2,l_axis2,p_axis,sign_phi)

xphi=xC+phi/pi*(xB-xC);
yphi=yC+phi/pi*(yB-yC);

R=[cos(sign_phi*phi+phi_0) -sin(sign_phi*phi+phi_0); ...
   sin(sign_phi*phi+phi_0)  cos(sign_phi*phi+phi_0)]';
D=[R -R*[xA;yA];...
   0 0 1];

Cp=[(R(2,:)*[xphi-xA;yphi-yA])/(R(1,:)*[xphi-xA;yphi-yA])^2 0 0;...
    0 0 -0.5;
    0 -0.5 0];

%The conic
C=D'*Cp*D;

%lBC=cross([xB;yB;1],[xC;yC;1]);
%d=getdirvec(C,xA,yA,xB,yB);
l_axis=cross([d;0],[xA;yA;1]); %X-axis of the boundary parabola
l0=cross(cross(l_axis,[0;0;1]),[X;Y;1]);

%Find the points on the line l0 that meet the conic
[x1,x2]=coniclinemeet(C,l0);
x1=real(x1)/real(x1(3));
x2=real(x2)/real(x2(3));

%Select the right solution
%
%The following criterion is quite complicated. It tests whether
%the direction has the correct sign AND the solution is in the region closed by the boundary parabolas
%Numerically this could perhaps still be improved: if x1 is close to the plane at
%infinity, the measurement of the sign of l_axis2'*x1 is questionable
TOL=10^(-6);
dB=getdirvec(C,xA,yA,xphi,yphi);
d1=getdirvec(C,xA,yA,x1(1),x1(2));
if sign(dB(1))==sign(d1(1)) & sign(dB(2))==sign(d1(2)) & ...
    ~(x1'*c1*x1>TOL & (sign(x1'*c1*x1)~=...
           sign(p_axis'*c1*p_axis)...
           & sign(l_axis2'*x1)==sign(l_axis2'*[xB;yB;1]))... 
           |...
           (x1'*c2*x1>TOL & sign(x1'*c2*x1)~=...
           sign(p_axis'*c2*p_axis)...
           & sign(l_axis2'*x1)==sign(l_axis2'*[xC;yC;1]))) 
%parab_y_axis=cross([[0 -1;1 0]*dB;0],[xA;yA;1]);
%if sign(parab_y_axis'*x1)==sign(parab_y_axis'*[xphi;yphi;1]);
%THIS DID NOT WORK AROUND THE P_AXIS SINCE THE TWO SOLUTIONS PROBABLY THEN HAVE
%NUMERICALLY THE SAME SIGN
x=x1;
else
    x=x2;
end

%the direction unit vector on the axis line
u=d; %[xB;yB]-[xC;yC];
%u=u/norm(u);

result=u'*[x(1)-X;x(2)-Y];
