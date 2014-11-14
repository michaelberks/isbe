function [ image_out ] = lloyd_max(image_in, plot_flag)
% lloyd_max Summary of this function goes here
%  Detailed explanation goes here

if nargin < 2
    plot_flag = 0;
end
[counts, x_range] = hist(double(image_in(:)), 256);
counts(end) = [];
x_range(end) = [];

goaty = max(counts);

L1 = min(x_range); L2 = max(x_range); y(1) = max(x_range);
i = 1;
display(['y1 = ',num2str(y(1))]);
while(true)
    i = i+1;
    y(i) = round((L1 + L2) / 2); %#ok
    display(['y',num2str(i),' = ',num2str(y(i))]);
    
    if plot_flag
        figure; bar(x_range, counts); hold on;
        plot([y(i) y(i)], [0 goaty], 'r');
        plot([L1 L1], [0 goaty], 'c');
        plot([L2 L2], [0 goaty], 'g');
    end
    
    x1 = x_range(x_range < y(i)); x2 = x_range(x_range >= y(i));
    Q1 = counts(x_range < y(i)); Q2 = counts(x_range >= y(i));
    
    L1 = round((Q1*x1')/sum(Q1)); L2 = round((Q2*x2')/sum(Q2));
    if y(i) == y(i-1), break
    end
    
end
display(['i = ', num2str(i)]);

image_out = image_in > y(end);

    
