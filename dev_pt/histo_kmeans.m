function [idx,c,sumD,D] = histo_kmeans(hists,k,distfcn,distmat)
% run K-means on histogram data with distance function specified

if nargin<3, distfcn = 'L1'; end
if nargin<4, distmat = []; end

f_debug = false;
if nargin==0 && nargout==0
	% debug mode
	f_debug = true;
	clc;
	
	hists = rand(3,200);
	hsum = sum(hists);
	hists = hists./hsum([1,1,1],:);
	k = 10;
end

% if histograms are in 3D format then reshape
if length(size(hists))==3
	hists = reshape(hists,[numel(hists(:,:,1)),size(hists,3)]);
end

N = size(hists,2);

% check that there are at least K examples
if N<k, error(sprintf('Need at least %i examples',k)); end

% choose K values at random
p = randperm(N);
c = hists(:,p(1:k));

idx0 = zeros(1,N);
sumD = zeros(1,k);
for it = 1:inf
	D = zeros(k,N);
	
	% compute distances between examples and centres
	for ic = 1:k
		for is = 1:N
			D(ic,is) = hist_dist(c(:,ic),hists(:,is),distfcn,distmat);
		end
	end
	
	% find nearest neighbours
	[minD,idx] = min(D);
	
	% replace centres with mean of histogram clusters and
	% compute sum of distances from each centre to its assigned samples
	for ic = 1:k
		c(:,ic) = mean(hists(:,idx==ic),2);
		sumD(ic) = sum(minD(idx==ic));
	end
	
	if f_debug
		% define basis for points and project
		theta = [0 1 2]*2*pi/3 + pi/2;
		basis = [cos(theta); sin(theta)];
		pts = basis*hists;
		c_pts = basis*c;

		figure(3); clf; hold on; box on;
			plot(basis(1,[1,2,3,1]),basis(2,[1,2,3,1]),'-','color',0.9*[1,1,1]);
			plot(pts(1,:),pts(2,:),'b.');
			plot(c_pts(1,:),c_pts(2,:),'ro');
			axis('equal',[-1,1,-0.6,1.1]);
		drawnow;
	end
	
	% if no change in assignments then terminate early
	if all(idx==idx0), break; end
	
	% otherwise, update idx0
	idx0 = idx;
end

if f_debug, clear; end

