function   ok = drawcircle(x,y,r,colour)

t=0:.01:2*pi;
line(x+r*sin(t),y+r*cos(t),'Color',colour,'LineWidth',2);  % r is the radius
ok=1;