function downsize_batch(file_names, folder_in, folder_out, height)

fid = fopen(file_names);
names = textscan(fid, '%s');
[x, y] = size(names{1});
for i=1:x,
    f_in = strcat(folder_in, names{1}{i});
    
    no_zero = 2 - floor(log10(i));
    zero_str = '';
    for j = 1:no_zero
        zero_str = strcat(zero_str, int2str(0));
    end
    
    f_out = strcat(folder_out, zero_str, int2str(i), '.bmp');
    
    i1 = imread(f_in);    
    i2 = 255 - i1; clear i1;
    [r, c, h] = size(i2);
    if r < c
        i2 = imrotate(i2, 90);
        temp = r; r = c; c = temp;
    end
    width = round(c*height/r);
    i3 = imresize(i2, [height, width], 'bilinear'); clear i2;
    imwrite(i3, f_out); clear i3;
    
end
fclose('all');