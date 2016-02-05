function show_ori_wheel(sz)

if nargin < 1
    sz = 100;
end
xy = repmat(-sz:sz, 2*sz+1, 1);
c = complex(-xy, xy');%.^2;
c((xy.^2 + xy'.^2) > sz^2) = 0; 
figure; imgray(complex2rgb(c));