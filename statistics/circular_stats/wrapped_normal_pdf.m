function [p] = wrapped_normal_pdf(theta, mu, rho)
% Defines the probability density function for a random variable theta that
% has a Wrapped Normal distribution WN(mu, rho);

%The exact PDF is computed as an infinite sum, however, for practical
%purposes,  we sum up to sum set of limits, [0 100] seems reasonable
% sigma2 = -2*log(rho);
% 
% 
% pdf_sum = 0;
% for k = -10:10
%     pdf_sum = pdf_sum + aux_sum(theta, mu, sigma2, k);
% end
% p = pdf_sum / sqrt(2*sigma2*pi);
% 
% function pdf_sum = aux_sum(theta, mu, sigma2, k)
%  
% pdf_sum = exp( -(theta - mu + 2*pi*k).^2 / (2*sigma2) );

pdf_sum = 0;
for k = 1:100
    pdf_sum = pdf_sum + aux_sum(theta, mu, rho, k);
end
p = (1 + 2*pdf_sum) / (2*pi);

function pdf_sum = aux_sum(theta, mu, rho, k)
 
pdf_sum = (rho^(k^2))*cos(k*(theta-mu));