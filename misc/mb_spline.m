function yi = mb_spline(x, y, xi, n)

if nargin < 4
    n = 3;
end

dims = size(y,2);
num_pts = length(xi);
yi = zeros(num_pts, dims);

for ii = 1:num_pts
    l = sum(x <= xi(ii));
    yi(ii,:) = compute_dik(x, y, xi(ii), l, n, n);
end

function dik = compute_dik(u, d, x, i, k, n)
    
    if k
        alpha = (x - u(i)) / (u(i+n+1-k) - u(i));
        dik = (1 - alpha)*compute_dik(u, d, x, i-1, k-1, n) + ...
            alpha*compute_dik(u, d, x, i, k-1, n);
    else
        dik = d(i);
    end