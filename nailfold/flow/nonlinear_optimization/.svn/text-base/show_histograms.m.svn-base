function show_histograms(p_test)

p_test.uu = zeros(size(p_test.uu));
p_test.vv = zeros(size(p_test.vv));

[v, sz_vec] = ncm_pack(p_test);
be = brightness_error_vec(v, sz_vec, p_test);

r = double(p_test.obsMask);
r(logical(r)) = be;

rmax = max(abs(r), [], 3);

button = 0;
while (button ~= 3)
    figure(9);
        imagesc(rmax);
        axis('image');
        
    [c,r,button] = ginput(1);
    z = r(round(cr(2)), round(cr(1)), :);
    z = z(z ~= 0);
    
    figure(8);
        hist(z, (-1:0.1:1));
end