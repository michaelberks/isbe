function [eout] = elle_eval_huber(X, alpha, beta, biv, bih)
% eout = elle_eval_huber(x,alpha,bet,biv,bih); */

r2 = sqrt(2);
eout = 0;

for i = 1:biv %(i=0; i<biv; i++) {
   for j = 1:bih %(j=0; j<bih; j++) {
       c = (j-1)*biv+i;
        if i<biv
            eout = eout + huber(X(c+1)-X(c), alpha);
            if j<bih
                eout = eout + huber((X(c+biv+1)-X(c))/r2, alpha); 
            end
        end
        if j<bih
            eout = eout + huber(X(c+biv)-X(c), alpha);
            if i>1  
                eout = eout + huber((X(c+biv-1)-X(c))/r2, alpha); 
            end
        end
   end
end
eout = eout * (beta/2); % now eout = nu*Huber(<x_neighbour_differences,alpha>); */


%-------------------------------------------------------------------------
function g = huber(x, alpha)
if (x<0) 
    x=-x;
end
if (x<alpha) 
    g = x.*x; 
else
    g = (2*x*alpha - alpha*alpha);
end

