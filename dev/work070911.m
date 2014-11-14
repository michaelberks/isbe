load C:\isbe\dev\files\u_files.mat
idx = randsample(101, 20);
r_files = u_files1(idx); 

[w1, o1] = optimise_weights(r_files, [0.5e-3 0.5e-3]);
save C:\isbe\dev\weights\opt1 w1 o1;
pack;

[w2, o2] = optimise_weights(r_files, [2.5e-3 0.5e-3]);
save C:\isbe\dev\weights\opt2 w2 o2;
pack;

[w3, o3] = optimise_weights(r_files, [0.5e-3 2.5e-3]);
save C:\isbe\dev\weights\opt3 w3 o3;
pack;