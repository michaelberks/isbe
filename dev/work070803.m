dim = 500;
%%
for ii = 1:N
    [dd Z t] = procrustes(reshape(m3, [], 2),...
        [X_shape(ii, 1:dim)' X_shape(ii, dim+1:end)']);
    X_pro4(ii, :) = [Z(:,1)', Z(:,2)'];
    %sq_diff(ii) = dd;
    %shape_scale(ii) = t.b;
    clear dd Z t;
end
    
%%
for ii = 1:N
    [dd Z t] = procrustes(ss,...
        [X_shape(ii, 1:dim)' X_shape(ii, dim+1:end)']);
    X_pro2(ii, :) = [Z(:,1)', Z(:,2)'];
    %sq_diff(ii) = dd;
    %shape_scale(ii) = t.b;
    clear dd Z t;
end
%%
for ii = 1:10
    figure('WindowStyle', 'docked');
    plot(X_pro1(ii,1:500), X_pro1(ii,501:end), 'b');
    hold on;
    plot(X_pro2(ii,1:500), X_pro2(ii,501:end), 'r:');
end
    