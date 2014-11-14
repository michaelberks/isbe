function flip_right(file_names, folder_in)

fid = fopen(file_names);
names = textscan(fid, '%s');
[x, y] = size(names{1});
for i=1:x,
    f_in = strcat(folder_in, names{1}{i});
    
    i1 = imread(f_in);    
    tform = maketform('affine',[-1 0 0; 0 1 0; 0 0 1]);
    i2 = imtransform(i1, tform); clear i1;
    imwrite(i2, f_in); clear i2;
    
end
fclose('all');