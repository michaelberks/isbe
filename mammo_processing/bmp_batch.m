function bmp_batch(file_names1, file_names2)

fid1 = fopen(file_names1);
names1 = textscan(fid1, '%s');

fid2 = fopen(file_names2);
names2 = textscan(fid2, '%s');

fclose('all');

[x, y] = size(names1{1});
for i=1:x,
    f_in = names1{1}{i};
    f_out = names2{1}{i};
    i1 = dicomread(f_in);
    i2 = im2uint8(i1); clear i1;
    imwrite(i2, f_out); clear i2;
end