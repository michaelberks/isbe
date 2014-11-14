function k = compute_kappa_statistic(oij)
% oi = sum(oij,1);
% oj = sum(oij,2);
% tc = sum(oi);
% 
% pr_e = (oi(1) / tc) * (oj(1) / tc) + (oi(2) / tc) * (oj(2) / tc);
% pr_a = (oij(1,1)+oij(2,2)) / tc;
% 
% k = (pr_a - pr_e) / (1 - pr_e);

a = oij(1,1);
b = oij(1,2);
c = oij(2,1);
d = oij(2,2);

k = 2*(a*d - b*c) / (2*a*d + b*b + c*c + a*(b+c) + d*(b+c));
