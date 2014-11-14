function intimg_diag = integral_image(img)
% function to compute the diagonal integral image (summed area table)
% as described in Lienhart and Maydt, ICIP 2002
%
% Note that the backward pass begins at row r-1 (not mentioned in the
% paper)

f_debug = false;
if nargin==0
	f_debug = true;
	img = rand(5,5);
end

% allocate space for integral image
[r,c] = size(img);
intimg_diag = zeros(r+1,c+2);

% compute integral image
intimg = zeros(r+1,c+2);
intimg(2:end,3:end) = cumsum(cumsum(img,1),2);

% forward pass
for i = 1:r
	intimg_diag(i+1,3:end) =	intimg(i+1,3:end) - ...
								intimg(i,3:end) + ...
								intimg_diag(i,2:end-1);
end

% backward pass
for i = r-1:-1:1
	intimg_diag(i+1,3:end) =	intimg_diag(i+1,3:end) + ...
								intimg_diag(i+2,2:end-1) - ...
								intimg_diag(i+1,1:end-2);

end

% trim extra rows and columns
intimg_diag = intimg_diag(2:end,3:end);

if f_debug
	clear;
end

