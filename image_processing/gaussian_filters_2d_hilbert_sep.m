function [filters] = gaussian_filters_2d_hilbert_sep(sigma, width)
% Auxiliary function to make the 2d gaussian derivative filters
if ~exist('width','var'), width = round(5*sigma); end

% Formula for the kernel taken from Appendix F of:
%   "The design and use of steerable filters",
%   W. T. Freeman and E. H. Adelson,
%   PAMI 13(9), September 1991

filters = zeros(4, 2*width+1);

sigmasq = sigma^2;
xx = (-width:width) / sqrt(2*sigmasq);
k = 1 / sqrt(2*pi*sigmasq);

filters(1, :) = k*sqrt(0.9780) *(-2.254*xx + xx.^3) .* exp(-xx.^2);
filters(2, :) = k*sqrt(0.9780) *(-0.7515 + xx.^2) .* exp(-xx.^2);
filters(3, :) = k*sqrt(0.9780) * xx .* exp(-xx.^2);
filters(4, :) = k*sqrt(0.9780) * exp(-xx.^2);

if nargout == 0
    
    H2a = filters(4, :)' * filters(1, :);
    H2b = filters(3, :)' * filters(2, :);
    H2c = filters(2, :)' * filters(3, :);
    H2d = filters(1, :)' * filters(4, :);
    
    figure;
    subplot(2,2,1); imgray(H2a); title('H_{2a}');
    subplot(2,2,2); imgray(H2b); title('H_{2b}');
    subplot(2,2,3); imgray(H2c); title('H_{2c}');
    subplot(2,2,4); imgray(H2d); title('H_{2d}');
end
    

