function [] = plot_cls_points(cls_results, line_style)

hold on;

for level = 3:3

    L = bwlabel(cls_results{level}.CLS, 8);

    for i = 1:max(L(:));
        [r c] = find(L == i);
        
        r = r*2^(level-1);
        c = c*2^(level-1);
        
        plot(c, r, line_style);        
    end
    
end