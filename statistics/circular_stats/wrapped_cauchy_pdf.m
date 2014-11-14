function [p] = wrapped_cauchy_pdf(theta, mu, rho)
% Defines the probability density function for a random variable theta that
% has a Wrapped cauchy distribution WN(mu, rho);

%The exact PDF is computed as an infinite sum, however, for practical
%purposes,  we sum up to sum set of limits, [-10 10] seems reasonable

pdf_sum = 0;
for k = 1:100
    pdf_sum = pdf_sum + aux_sum(theta, mu, rho, k);
end
p = (1 + 2*pdf_sum) / (2*pi);

function pdf_sum = aux_sum(theta, mu, rho, k)
 
pdf_sum = (rho^k)*cos(k*(theta-mu));