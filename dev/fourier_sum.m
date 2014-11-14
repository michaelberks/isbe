function s_x = fourier_sum(x, p, pheta, lim, if_plot)



% Fourier series sum
if if_plot; figure; hold on; end

s_x = zeros(1,length(x));
for n = 0:lim
    
    sum_term = sin( (2*n + 1).*x + pheta ) / (2*n + 1)^p;
    s_x = s_x + sum_term;
    if if_plot; plot(x, sum_term, 'r:'); end
    
end

if if_plot; plot(x, s_x, 'b-'); end
    
    