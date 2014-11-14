function [output1, output2] = localimagemaximanocomments(input1, input2, input3, input4)
%LOCAL_IMAGE_MAXIMA_NO_COMMENTS does stuff
[r c] = size(input1);
if nargin < 2,input2 = 1;end
if nargin < 3 || isempty(input3),input3 = true(r, c);end
if nargin < 4,input4 = -inf;
end
input3 = input3 & (input1 > input4);
input1 = padarray(input1, [1 1], -inf);
A = input3;
for i = -1:1
for j = -1:1
if ~i && ~j
continue
end
A = A & (input1(2:end-1, 2:end-1) >= input1((2:end-1)+i, (2:end-1)+j));
end
end
input1 = input1(2:end-1, 2:end-1);
B = find(A);
[a b] = ind2sub([r c], B);
[sm si] = sort(input1(B), 'descend');
b2 = b(si);
a2 = a(si);
n = length(sm);
k = true(n, 1);
for i = 1:n
if k(i)
x = b2(i);y = a2(i);
k = k & ((b2-x).^2 + (a2-y).^2 > input2^2);k(i) = true;
end
end
output1 = [b2(k) a2(k)];
if nargout > 1
output2 = input1(B(si(k)));
end
