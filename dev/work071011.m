figure; hold on;
for ii = 1:101
    plot(r_pts(1,:,ii), r_pts(2,:,ii));
    plot(r_pts(1,1,ii), r_pts(2,1,ii), 'r.', 'MarkerSize', 10);
end
%%
figure; hold on;
for ii = 1:101
    plot(shapes1(1,:,ii), shapes1(2,:,ii));
    plot(shapes1(1,1,ii), shapes1(2,1,ii), 'r.', 'MarkerSize', 10);
end
figure; hold on;
for ii = 1:101
    plot(shapes2(1,:,ii), shapes2(2,:,ii));
    plot(shapes2(1,1,ii), shapes2(2,1,ii), 'r.', 'MarkerSize', 10);
end
%%
figure; hold on;
for ii = 1:101
    plot(shapes_p(ii, 1:end/2), shapes_p(ii, end/2+1:end));
    plot(shapes_p(ii, 1), shapes_p(ii, end/2+1), 'r.', 'MarkerSize', 10);
end
figure; hold on;
for ii = 1:101
    plot(shapes_pro(ii, 1:end/2), shapes_pro(ii, end/2+1:end));
    plot(shapes_pro(ii, 1), shapes_pro(ii, end/2+1), 'r.', 'MarkerSize', 10);
end
%%
figure; hold on;
for ii = 1:101
    plot(a_shapes(ii, 1:end/2), a_shapes(ii, end/2+1:end));
    plot(a_shapes(ii, 1), a_shapes(ii, end/2+1), 'r.', 'MarkerSize', 10);
end