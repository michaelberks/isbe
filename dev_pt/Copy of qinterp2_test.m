clc;

m = 10; n = 12;
img = repmat(linspace(0,255,n),[m,1]) + ...
			repmat(linspace(0,255,m)',[1,n]);
img = img/max(img(:));

% sample points as a grid
[x,y] = meshgrid(1:0.25:n,1:0.25:m);

interptype = 'linear';
n_tests = 1000;

tic; 
	for i = 1:n_tests
		img3 = interp2(img,x,y,'linear'); 
	end
toc;

switch interptype
	case 'linear',
		tic;
			for i = 1:n_tests
				% get integer parts
				xf = floor(x); xc = ceil(x);
				yf = floor(y); yc = ceil(y);

				% get fractional parts
				xr = x-xf; yr = y-yf;

				q1	= (1-xr(:)) .* (1-yr(:)) .* img(1+(xf(:)-1)*m+(yf(:)-1));
				q2	=    xr(:)  .* (1-yr(:)) .* img(1+(xc(:)-1)*m+(yf(:)-1));
				q3	= (1-xr(:)) .*    yr(:)  .* img(1+(xf(:)-1)*m+(yc(:)-1));
				q4	=    xr(:)  .*    yr(:)  .* img(1+(xc(:)-1)*m+(yc(:)-1));
				img2 = reshape(q1+q2+q3+q4,size(x));
			end
		toc;
		
	case 'cubic',

end

figure(1); 
	subplot(3,1,1); imshow(img);
	subplot(3,1,2); imshow(img2);
	subplot(3,1,3); imshow(img3);
