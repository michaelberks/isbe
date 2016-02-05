
function [result,dresultda]=intbody(x,a)

%Result of the integral sqrt(ax^2+bx+c)
%result=(b+2*a*x)/(4*a)*sqrt(a*x^2+b*x+c)+(4*a*c-b^2)/(8*a^(3/2)) * ...
%    log(abs(2*a*x+b+2*sqrt(a*(a*x^2+b*x+c))));

%Result of the integral sqrt(a^2+x^2)
result=1/2*(x*sqrt(a^2+x^2)+a^2*log(x+sqrt(a^2+x^2)));
dresultda=1/2*(x*a*(a^2+x^2)^(-1/2)+2*a*log(x+sqrt(a^2+x^2))+a^3*(x+sqrt(a^2+x^2))^(-1)*(a^2+x^2)^(-1/2));
