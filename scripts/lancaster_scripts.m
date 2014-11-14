
tree = imresize(imread('K:\isbe\conferences_and_symposia\lancaster_wavelets_jan2012\figures\tree_gray.bmp'), 0.5, 'bilinear');
for sigma = [1 2 4 8 16 32]
    
%     [line_strength, orientation] = gaussian_1st_derivative_gradient2(tree, sigma);
% 
%     imwrite(complex2rgb(line_strength .* exp(2*1i*orientation)),...
%         ['K:\isbe\conferences_and_symposia\lancaster_wavelets_jan2012\figures\tree_g1d_' num2str(sigma) '.bmp']);
    
    figure; image(complex2rgb(line_strength .* exp(2*1i*orientation))); axis image;
    [line_strength, orientation] = gaussian_2nd_derivative_line(tree, sigma);

    figure; image(complex2rgb(line_strength .* exp(2*1i*orientation))); axis image;
    imwrite(complex2rgb(line_strength .* exp(2*1i*orientation)),...
        ['K:\isbe\conferences_and_symposia\lancaster_wavelets_jan2012\figures\tree_g2d_' num2str(sigma) '.bmp']);
end
%%
%%
%-------------------------------------------------------------------------
nodes = [0 0];
next_nodes = [];
len = 10;
level = 0;
figure; hold on;
%while ~isempty(nodes);
for level = 0:4
    level_x = 2^(6 - level);
    level_y = len;
    level_p = 1 - 2^(-level/4);
    for n = 1:size(nodes,1);
        plot(nodes(n,1), nodes(n,2), 'o', 'markersize', 20, 'markeredgecolor', 'r', 'markerfacecolor', 'c');
        if rand > level_p
            child1 = [nodes(n,1) + level_x nodes(n,2) + level_y];
        
            plot([nodes(n,1) child1(1)], [nodes(n,2) child1(2)]);
            next_nodes = [next_nodes; child1];
            plot(child1(1), child1(2), 'o', 'markersize', 20, 'markeredgecolor', 'r', 'markerfacecolor', 'c');
        end
        if rand > level_p
            child2 = [nodes(n,1) - level_x nodes(n,2) + level_y];    
            plot([nodes(n,1) child2(1)], [nodes(n,2) child2(2)]);
            next_nodes = [next_nodes; child2];
            plot(child2(1), child2(2), 'o', 'markersize', 20, 'markeredgecolor', 'r', 'markerfacecolor', 'c');
        end
    end
    nodes = next_nodes;
    next_nodes = [];
    %level = level + 1;
end