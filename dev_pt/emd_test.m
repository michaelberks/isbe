clc;

% number of bins in each
m = 3; n = 3;
n = m;

% generate two histograms
p = rand(1,m); p = p/sum(p(:));
q = rand(1,n); q = q/sum(q(:));
% p = ones(1,m); p = p/sum(p(:));
% q = ones(1,n); q = q/sum(q(:));


% determine distances between cells
cp = (1:m)/m - 1/(2*m);
cq = (1:n)/n - 1/(2*n);
[xx,yy] = meshgrid(cq,cp);
% alpha = 1; d_ij = 1 - exp(-alpha*abs(xx-yy)); % metric (needs improving)
d_ij = abs(xx*n-yy*m);

c_i = zeros(m,n,m); for i = 1:m, c_i(i,:,i) = 1; end
c_j = zeros(m,n,n); for j = 1:n, c_j(:,j,j) = 1; end
c = ones(m,n);

% set up linear programming problem
A	= []; b = [];
% A	= [	A; eye(m*n,m*n) ]; b = [ b; zeros(m*n,1) ];
A	= [ A; reshape(c_i,[m*n,m])' ]; b = [ b; p(:) ];
A	= [ A; reshape(c_j,[m*n,n])' ]; b = [ b; q(:) ];
Aeq	= []; beq = [];
Aeq	= [ Aeq; ones(1,m*n) ]; beq = [ beq; min(sum(p(:)),sum(q(:))) ];
LB	= zeros(m*n,1);
UB	= [];
X0	= zeros(m*n,1);

tic;
f_emd = linprog(d_ij(:),A,b,Aeq,beq,LB,UB);
d_emd = f_emd'*d_ij(:)/sum(f_emd(:));
toc;

d2 = hist_dist(p,q,'emd');

if m==n
	tic;
	d_match = hist_dist(p,q,'match');
	toc;
end

display([p(:) reshape(f_emd,[m,n]); nan q]);
display([d_emd d2 d_match]);


