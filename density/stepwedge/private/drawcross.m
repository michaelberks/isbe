function   ok = drawcross(x,y,size,colour)
% draws a cross
% X is vector containing x coords of opposite corners 
% Y is vector containing y coords of opposite corners

line([x-size x+size],[y-size y+size],'Color',colour,'LineWidth',2);
line([x+size x-size],[y-size y+size],'Color',colour,'LineWidth',2);
   
    
ok=1;