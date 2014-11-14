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
		img2 = interp2(img,x,y,'linear'); 
	end
toc;

profile clear; profile on;
tic; 
	for i = 1:n_tests
		img3 = qinterp2(img,x,y,'linear'); 
	end
toc;
profile report

figure(1); 
	subplot(3,1,1); imshow(img);
	subplot(3,1,2); imshow(img2);
	subplot(3,1,3); imshow(img3);
	
err = mean(abs(img2(:)-img3(:)))
