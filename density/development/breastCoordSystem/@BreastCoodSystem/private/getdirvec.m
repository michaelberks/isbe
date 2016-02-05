function d_A=getdirvec(C,xA,yA,xB,yB)

%Computes the direction vector of the parabola C at (xA,yA) to the
%direction of another point (xB,yB) on the parabola

zA=x2z(xA,yA,C);
zB=x2z(xB,yB,C);

[tmp1,tmp2,dx,dy]=z2x(zA,C);

d_A=[dx;dy];
d_A=sign(zB-zA)*d_A/norm(d_A);

