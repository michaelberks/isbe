function d = hist_dist(h1,h2,distname,varargin)

if nargin==0 && nargout==0
	% debug mode
	clc;
	n = 10; sig = 1;
	h1 = rand(1,n); h2 = abs(h1 + sig*randn(size(h1)));
	distnames = {...
		'bhattacharyya','l1','l2','lmax',...
		'kl','jeffrey',...
		'intersection','chisquared',...
		'quadratic','kolmogorov-smirnov','match','emd'};
	for i = 1:length(distnames)
		fprintf('%20s(h1,h2) = %f\n',...
				distnames{i},hist_dist(h1,h2,distnames{i}));
	end
end

% default inputs
if nargin<3, distname = 'bhattacharyya'; end

% default outputs
d = NaN;

% check dimensions match
n = numel(h1);
if (numel(h2) ~= n) && ....
	~any(strcmpi(distname,{'emd','earthmovers'}))
	
	error('Histograms must have same dimensionality');
end

% check for empty histograms
if all(h1(:)==0) || all(h2(:)==0), return; end

% check for histograms with negative entries
if any(h1(:)<0) || any(h2(:)<0)
	error('Histograms must be positive for all i');
end


% check applicability for certain methods
switch lower(distname)
	case {'match','kolmogorov-smirnov'}
		if (sum(size(h1)>1) > 1)
			error([distname,' distance measure not defined for >1 dimension']);
		end
end


% normalize
switch lower(distname)
	case 'bhattacharyya',
		% normalize to norm=1
		h1 = h1(:); h2 = h2(:);
		n1 = norm(h1(:)); if abs(n1-1)>1e-6, h1 = h1(:)/n1; end
		n2 = norm(h2(:)); if abs(n2-1)>1e-6, h2 = h2(:)/n2; end
	
	case {'l1','manhattan','cityblock',...
		  'l2','euclidean',...
		  'lp','minkowski',...
		  'lmax','linf','chebyshev',...
		  'kl','jeffrey','chisquared',...
		  'quadratic',...
		  'match','kolmogorov-smirnov'}
		% normalize to sum=1
		h1 = h1(:); h2 = h2(:);
		s1 = sum(h1(:)); if abs(s1-1)>1e-6, h1 = h1(:)/s1(:); end
		s2 = sum(h2(:)); if abs(s2-1)>1e-6, h2 = h2(:)/s2(:); end
		
	case {'emdn'}
		% normalize to sum=1
		h1 = h1(:); h2 = h2(:);
		s1 = sum(h1(:)); if abs(s1-1)>1e-6, h1 = h1(:)/s1(:); end
		s2 = sum(h2(:)); if abs(s2-1)>1e-6, h2 = h2(:)/s2(:); end
		distname = 'emd';
end


% compute normalized distance (not similarity)
switch lower(distname)
	case 'bhattacharyya',
		d = 1 - h1'*h2;
		
	case {'l1','manhattan','cityblock'}
		% max L1 dist = 2^(1/1) = 2
		diff = h1-h2;
		d = sum(abs(diff)) / 2;
		
	case {'l2','euclidean'}
		% max L2 dist = 2^(1/2) = sqrt(2)
		diff = h1-h2;
		d = sqrt(sum(diff.*diff)) / sqrt(2);
		
	case {'lp','minkowski'}
		% max Lp dist = 2^(1/p)
		if length(varargin)<1
			error('Lp distance requires a value for p');
		end
		p = varargin{1};
		diff = h1-h2;
		if isinf(p),	d = max(abs(diff));
		else			d = sum(abs(diff).^p).^(1/p) / (2^(1/p));
		end
		
	case {'lmax','linf','chebyshev'}
		% max Linf dist = 2^(1/inf) = 2^0 = 1
		diff = h1-h2;
		d = max(abs(diff));
		
	case {'kl'}
		inds = (h1>0);
		if any(h2(inds)==0), return; end % not defined - return NaN
		d = sum( h1(inds) .* log(h1(inds)./h2(inds)) );

	case {'jeffrey'}
		m = (h1+h2)/2;
		inds = (h1>0) & (h2>0);
		if any(m(inds)==0), return; end % not defined - return NaN
		d = sum( h1(inds) .* log(h1(inds)./m(inds)) + ...
				 h2(inds) .* log(h2(inds)./m(inds)) );
			 
	case {'chisquared'}
		m = (h1+h2)/2;
		diff = h1(:)-m(:);
		d = sum( diff.*diff ./ m(:) );
		
	case {'intersection'}
		d = 1 - sum(min(h1(:),h2(:)))./sum(h2(:));

	case {'quadratic'}
		if ~isempty(varargin)
			d_ij = varargin{1};
		else
			warning('Using 1-d_ij/d_max for quadratic distance');
			[xx,yy] = meshgrid(1:n,1:n);
			d_ij = abs(xx-yy);
		end
		
		d_max = max(d_ij(:));
		A = 1 - d_ij/d_max;
		
		diff = h1-h2;
		d = sqrt(diff'*A*diff);
		
	case {'match'}
		c1 = cumsum(h1); c2 = cumsum(h2);
		diff = c1-c2;
		d = sum(abs(diff));
		
	case {'kolmogorov-smirnov'}
		c1 = cumsum(h1); c2 = cumsum(h2);
		diff = c1-c2;
		d = max(abs(diff));
		
	case {'emd','earthmovers'}
		% this is equivalent to 'match' if both histograms are the same 
		% size and d_ij=abs(pi-qi) and sum(h1)=sum(h2)=1
		h1 = h1(:); h2 = h2(:);
		m = length(h1); n = length(h2);
		
		if ~isempty(varargin)
			d_ij = varargin{1};
		else
			% just one choice of metric (though not an especially good one)
			c1 = (1:m)/m - 1/(2*m);
			c2 = (1:n)/n - 1/(2*n);
			[xx,yy] = meshgrid(c2,c1);
% 			alpha = 1; d_ij = 1 - exp(-alpha*abs(xx-yy));
			d_ij = abs(xx*n-yy*m);
		end

		% matrices representing constraints
		c_i = zeros(m,n,m); for i = 1:m, c_i(i,:,i) = 1; end
		c_j = zeros(m,n,n); for j = 1:n, c_j(:,j,j) = 1; end

		% set up linear programming problem
		A	= []; b = [];
		A	= [ A; reshape(c_i,[m*n,m])' ]; b = [ b; h1(:) ];
		A	= [ A; reshape(c_j,[m*n,n])' ]; b = [ b; h2(:) ];
		Aeq	= []; beq = [];
		Aeq	= [ Aeq; ones(1,m*n) ]; beq = [ beq; min(sum(h1(:)),sum(h2(:))) ];
		LB	= zeros(m*n,1);
		UB	= [];

		% solve linear programming problem and compute distance
		% note that this required the optimization toolbox
		opts = optimset('linprog');
			opts = optimset(opts,'display','none');
		f_emd = linprog(d_ij(:),A,b,Aeq,beq,LB,UB,[],opts);
		d = f_emd'*d_ij(:)/sum(f_emd(:));
		
	otherwise,
		error(['Unknown distance measure: ',distname]);
end

if nargin==0
	clear;
end
