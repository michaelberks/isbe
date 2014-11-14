function no_black_batch(file_names_o, file_names_f, safe)

if nargin < 3
    safe = 1;
end

fid_o = fopen(file_names_o);
names_o = textscan(fid_o, '%s');

fid_f = fopen(file_names_f);
names_f = textscan(fid_f, '%s');

fclose('all');

[x, y] = size(names_o{1});
for i=1:x,
    f_orig = names_o{1}{i};
    f_filt = names_f{1}{i};
    no_black(f_orig, f_filt, safe);
end