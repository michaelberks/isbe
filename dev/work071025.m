%%
%
% Comute the difference between means of leave-one-out models
means_diff = zeros(101,1000);
for ii = 1:101
    load(['C:\isbe\dev\models\loo_u1\model', zerostr(ii,3)]);
    means_diff(ii,:) = mass_model.mean_shape;
    clear mass_model;
end
%%
figure; axis equal; hold on;
for ii = 1:101
    plot(means_diff(ii,1:500), means_diff(ii,501:end));
end
%%
% Comute the difference between texture vectors of leave-one-out models
for ii = 2:101
    load(['C:\isbe\dev\models\loo_u1\model', zerostr(ii,3)]);
    tex_diff1(ii,:) = mass_model.X_tex(1,:);
    clear mass_model;
end
for ii = 1:100
    load(['C:\isbe\dev\models\loo_u1\model', zerostr(ii,3)]);
    tex_diff2(ii,:) = mass_model.X_tex(end,:);
    clear mass_model;
end
    