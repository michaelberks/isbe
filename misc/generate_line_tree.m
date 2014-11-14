function line_mask = generate_line_tree(breast_mask, start_pt, theta, line_mask, line_sum)

len = rand*100;
theta = theta + rand*pi/6 - pi/12;

end_pt(1) = round(start_pt(1) - len * cos(theta));
end_pt(2) = round(start_pt(2) - len * sin(theta));

    
if  all(end_pt >= 1) && end_pt(2) <= size(breast_mask, 1) && breast_mask(end_pt(2), end_pt(1))
    
    %plot([start_pt(1) end_pt(1)], [start_pt(2) end_pt(2)], color);
    
    if start_pt(1) > end_pt(1)
        cols = start_pt(1) - end_pt(1) + 1;
        sx = cols;
        sc = end_pt(1);
        
    else
        cols = end_pt(1) - start_pt(1) + 1;
        sx = 1;
        sc = start_pt(1);
    end
    
    if start_pt(2) > end_pt(2)
        rows = start_pt(2) - end_pt(2) + 1;
        sy = rows;
        sr = end_pt(2);
        
    else
        rows = end_pt(2) - start_pt(2) + 1;
        sy = 1;
        sr = start_pt(2);
    end
        
    [line_patch] = create_rect_bar(1, 1, 180*(pi-theta)/pi, rows, cols, sx, sy);
    line_mask(sr:sr+rows-1, sc:sc+cols-1) = line_mask(sr:sr+rows-1, sc:sc+cols-1) | ...
        (line_patch & breast_mask(sr:sr+rows-1, sc:sc+cols-1));
    
    if sum(line_mask(:)) > line_sum
        return;
    else
        if rand > .5
            line_mask = generate_line_tree(breast_mask, end_pt, theta, line_mask, line_sum);
        end
        if rand > .5
            line_mask = generate_line_tree(breast_mask, end_pt, theta, line_mask, line_sum);
        end
    end
    
end
        
