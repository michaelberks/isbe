function img = kars_test_image

% parameters
s = 1;			% relative size of image with respect to circle
r1 = 0.6;		% inner radius
r2 = 0.9;		% outer radius
N = s*128;		% size of the image
n_lines = 12;	% number of lines

% create grid of radius and theta
img = zeros(N,N);
[xx,yy] = meshgrid(linspace(-s,s,N),linspace(-s,s,N));
radius = sqrt(xx.*xx+yy.*yy);
theta = mod(atan2(yy,xx),pi);

% generate image with n_lines lines per half circle
trng = linspace(0,pi,n_lines+1); 
trng = (trng(1:end-1)+trng(2:end))/2;
for t = trng
	d = theta-t;
	[ignore,indsx] = min(abs(d),[],1);
	inds = sub2ind([N,N],indsx,1:N);
	img(inds) = t;
	[ignore,indsy] = min(abs(d),[],2);
	inds = sub2ind([N,N],1:N,indsy');
	img(inds) = t;
end

% blank out anything outside of the 'doughnut'
img(radius<r1) = 0;
img(radius>r2) = 0;

