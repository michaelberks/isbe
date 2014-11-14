function p = merge_structs(p1, p2)
% Merge structures p1 and p2, with p2 overwriting any existing fields in p1.

p = p1;

fnames = fieldnames(p2);
for i = 1:length(fnames)
    p.(fnames{i}) = p2.(fnames{i});
end
