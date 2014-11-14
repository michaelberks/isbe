nailfold = u_load('C:\isbe\nailfold\data\rsa_study\test\images\11200c.mat');
figure; imgray(nailfold);
circles_patch = nailfold(700:900, 600:900);
clear nailfold;

nailfold = u_load('C:\isbe\nailfold\data\rsa_study\test\images\10598c.mat');
figure; imgray(nailfold);
vessels_patch = nailfold(300:500, 800:1100);

%%
figure; imgray(circles_patch);
[discard_map] = remove_circles(circles_patch, [2 4], 6, 2, pi/10);
[cy cx] = find(discard_map);
plot(cx, cy, 'r.');

figure; imgray(vessels_patch);
[discard_map] = remove_circles(vessels_patch, [2 4], 6, 2, pi/10);
[cy cx] = find(discard_map);
plot(cx, cy, 'r.');
%%
im_list = dir('C:\isbe\nailfold\data\rsa_study\test\images\*.mat');
%%
for i_im = 7:20
    nailfold = u_load(['C:\isbe\nailfold\data\rsa_study\test\images\' im_list(i_im).name]);
    nailfold = imresize(nailfold, 0.5);
    figure; imgray(nailfold);
    [discard_map] = remove_circles(nailfold, [1 2], 6, 2, pi/10);
    [cy cx] = find(discard_map);
    plot(cx, cy, 'r.');
end
%%
im_list = dir('C:\isbe\nailfold\data\rsa_study\test\predictions\detection\rf_classification\257273\*.mat');
for i_im = [1:5 7:20]
    nailfold = u_load(['C:\isbe\nailfold\data\rsa_study\test\predictions\detection\rf_classification\257273\' im_list(i_im).name]);
    nailfold = imresize(nailfold, 0.5);
    figure; imgray(nailfold);
    [discard_map] = remove_circles(nailfold, [1 2], 8, 0.05, pi/12);
    [cy cx] = find(discard_map);
    plot(cx, cy, 'r.');
end