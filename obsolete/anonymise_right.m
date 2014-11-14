function anonymise_right(file_names, folder_in)

fid = fopen(file_names);
names = textscan(fid, '%s');
[x, y] = size(names{1});
for i=1:x,
    f_in = strcat(folder_in, names{1}{i});
    
    i1 = imread(f_in);
    grey_val = mean(mean(i1(1:350, end-85:end)));
    i1(1:350, end-85:end) = grey_val;
    imwrite(i1, f_in); clear i1;
    
end
fclose('all');