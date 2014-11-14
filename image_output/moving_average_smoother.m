function [sorted_x, moving_avg] = moving_average_smoother(x, y, width, N)

x = x(:);
y = y(:);

n_pts = length(x);
if nargin < 3
    width = round(n_pts / 100);
end
if ~rem(width,2)
    width = width + 1;
end
hw = (width-1)/2;

%Sort y with respect to sorted x
[sorted_x, inds] = sort(x);
sorted_y = y(inds);

%take moving average of y
smoother = gaussian_filters_1d(hw/5, hw);
smoother = smoother ./ sum(smoother);
moving_avg = conv(sorted_y, smoother, 'valid');

%discard start and end of x

sorted_x([1:hw end-hw+1:end]) = [];
n_pts = n_pts - width + 1;

if nargin > 3 && N < n_pts
    x = 1:n_pts;
    xi = linspace(1, n_pts, N);
    sorted_x = interp1(x, sorted_x, xi, 'linear');
    moving_avg = interp1(x, moving_avg, xi, 'linear');
end
