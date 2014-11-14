cur_pts = my_snake;
alpha = 1;
beta = 1;
feat_img = normal_med;

%alpha sum is total distance
alpha_sum = alpha*sum(sum(diff(cur_pts).^2, 2))
beta_sum = beta*sum(sum(diff(cur_pts, 2).^2, 2))
f_sum = sum(feat_img(sub2ind(size(feat_img), cur_pts(:,2), cur_pts(:,1))))

total_sum = alpha_sum + beta_sum - f_sum

%%

%%
normal_med = medfilt2(round(normal_p'));
norm_edge = imfilter(normal_med, [1 2 4 0 -4 -2 -1]', 'symmetric');
nms = (norm_edge(2:end-1,:) >= norm_edge(1:end-2,:)) & (norm_edge(2:end-1,:) >= norm_edge(3:end,:));
nms = [zeros(1,200); nms; zeros(1,200)];
norm_nms = norm_edge;
norm_nms(~nms) = 0;
norm_nms(norm_nms < 0) = 0;
%%
len = 200;
normal_med = medfilt2(normal_p);
norm_edge = imfilter(normal_med, [1 2 4 0 -4 -2 -1], 'replicate');

% nms = (norm_edge(:,2:end-1) < norm_edge(:,1:end-2)) | (norm_edge(:,2:end-1) < norm_edge(:,3:end));
% nms = [true(len,1) nms true(len,1)];
% norm_nms = norm_edge;

nms = (norm_edge(:,2:end-1) >= norm_edge(:,1:end-2)) & (norm_edge(:,2:end-1) >= norm_edge(:,3:end));
nms = [zeros(len,1) nms zeros(len,1)];
norm_nms = norm_edge;

norm_nms(~nms) = 0;
norm_nms(norm_nms < 0) = 0;
norm_nms = norm_nms';
%%
initial_line = [(1:200)' 25*ones(200,1)];
go_on = true;
ii = 2;
[snake_pnts,e] = snake(initial_line, 10, 10, 0, 0, 20, 1, norm_nms);
while go_on
    [snake_pnts,e(ii)] = snake(snake_pnts, 10, 10, 0, 0, 20, 1, norm_nms);
    go_on = e(ii) < e(ii-1);
    ii = ii+1;
end


