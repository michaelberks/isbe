clc;

figpath = 'u:/matlab/figs/mammography/';

clear forest;

% load([asymmetryroot,'data\line_orientation_rfs/pc20110211T143503/random_forest01.mat']);
% forest = random_forest;

if ~exist('forest','var')
	forest_root		= [asymmetryroot,'data/line_orientation_rfs/'];
	
	forest_dir		= dir([forest_root,'pc*']);
	forest_job		= [forest_dir(end).name,'/'];
	
	forest_dir		= dir([forest_root,forest_job,'/random_forest*.mat']);
	forest_fname	= [forest_root,forest_job,forest_dir(end).name];
	forest				= u_load(forest_fname);
	
	% load training data
	load([forest_root,forest_job,'01_trees/traindata.mat']);
end

in		= X(:,1:30).*exp(complex(0,X(:,31:end)));
in		= reshape(in,[size(in,1),6,5]);
out		= angle(y);

% set up labels for variables
labels = {};
for level = 1:5
	for subband = 1:6
		labels{subband,level} = sprintf('L%i:S%i',level,subband);
	end
end

ncols	= 256; t = linspace(0,2*pi,ncols);
cmap	= [(cos(t)+1)/2; (sin(t)+1)/2; zeros(1,ncols)]';

% default parameters
axlim = [];

%% get inputs in pairs
pair = 'ph_mag';
switch lower(pair)
	case 're_im'
		% real vs imag
		in = in(:,:); labels = labels(1:end);
		in = cat(3,real(in),imag(in)); 
		labels = cat(1,labels,labels);
		prefix = cat(1,repmat({'Re'},[1,size(in,2)]),repmat({'Im'},[1,size(in,2)]));
		fname = 'dt_scatter_Re_Im';
		sbsz = [5,6];
		
	case 'ph_mag',
		% phase vs magnitude
% 		in = in(:,:); labels = labels(1:end);
% 		in = cat(3,angle(in),log(abs(in)));
		in = in(:,:); labels = labels(1:end);
		in1 = angle(in); in2 = log(abs(in));
		labels = cat(1,labels,labels);
		prefix = cat(1,repmat({'Ph'},[1,size(in,2)]),repmat({'Mag'},[1,size(in,2)]));
		fname = 'dt_scatter_Ph_Mag';
		sbsz = [5,6];
		axlim = [-pi,pi,-3,3];

	case 'mag_mag',
		% mag vs mag, phase vs phase
		in = in(:,[1,6,2,5,3,4],:); in = in(:,:); % put similar orientations next to each other
		labels = labels([1,6,2,5,3,4],:); labels = labels(1:end);
		in = cat(3,[abs(in(:,1:2:end)) angle(in(:,1:2:end))],...
							 [abs(in(:,2:2:end)) angle(in(:,2:2:end))]); 
		in = reshape(in,[size(in,1)*3,size(in,2)/3,2]);
		in = [in(:,1:size(in,2)/2,:); in(:,1+size(in,2)/2:end,:)];
		in = reshape(in,[size(in,1)/6,size(in,2)*6,2]);
		labels = cat(2,reshape(labels,[2,length(labels)/2]),reshape(labels,[2,length(labels)/2]));
		prefix = cat(1,repmat({'Mag'},[6,size(in,2)/6]),repmat({'Ph'},[6,size(in,2)/6]));
		prefix = reshape(prefix,[2,size(in,2)]);
		fname = 'dt_scatter_Mag_Mag';
		sbsz = [5,6];

	case 'mag_mag2',
		% mag vs mag, phase vs phase
		in1 = [abs(in(:,:,2:end)) 2*angle(in(:,:,2:end))]; labels1 = labels(:,2:end);
		in2 = [abs(in(:,:,1:end-1)) angle(in(:,:,1:end-1))]; labels2 = labels(:,1:end-1);
		
		sbsz = [8,6];
	
		fname = 'dt_scatter_Mag_Mag2';
		
	case 'relmag',
		% mag vs mag, phase vs phase
		in1 = [log(abs(in(:,[2,3,4,5,6,1],:))./abs(in(:,[1,2,3,4,5,6],:)))]; labels1 = labels(:,2:end);
		in2 = [angle(in)]; labels2 = labels(:,1:end-1);
		
		sbsz = [5,6];
		axlim = [-5,5,-pi,pi];
	
		fname = 'dt_scatter_RelMag';
end

figure(1); clf; colormap(cmap);
d = 1; inds = reshape(1:prod(sbsz),sbsz(end:-1:1))';
for x1 = 1:sbsz(1)
	for x2 = 1:sbsz(2)
		mysubplot(sbsz,inds(x1,x2)); hold on;
		scatter(in1(:,d), in2(:,d), 4, out, '.');
		if ~isempty(axlim), axis(axlim); end
% 		axis('equal'); 
		set(gca,'box','on');
% 		plot(1e6*[-1,1],1e6*[-1,1],'-','color',0.7*[1,1,1]);
% 		xlabel(sprintf('%s(%s)',prefix{1,d},labels{1,d}));
% 		ylabel(sprintf('%s(%s)',prefix{2,d},labels{2,d}));
		d = d+1;
	end
end

set(1,'paperposition',[0,0,40,30]);
exportfig([figpath,fname],{'fig','eps','pdf','png'});

return


D = forest.D;
varhst = zeros(100,D);

vardisp = zeros(100,D);
vardispcount = zeros(100,D);

parmat = zeros(D,D);

width = [];
contrast = [];
orientation = [];
squash = [];

max_level = 0;

% figure(1); clf; hold on;
for t = 1:length(forest.trees)
	% load parameters
	load([forest.tree_root,forest.tree_dir,...
				sprintf('../line_parameters/01/parameters%03d.mat',t)]);
	width = [width [parameters.width]];
	contrast = [contrast [parameters.contrast]];
	orientation = [orientation [parameters.orientation]];
	squash = [squash [parameters.squash]];
	
			
	% load tree
	load([forest.tree_root,forest.tree_dir,forest.trees{t}]);
	leaves = find(all(tree.children==0,2));

	% sort leaves in order of output
	[sorted,inds] = sort(angle(tree.class(leaves))/2);
	leaves = leaves(inds);

% 	figure(1);
% 		for L = 1:length(leaves)
% 			leaf = leaves(L);
% 			n = length(tree.outputs{leaf});
% 			plot(L*ones(1,n),angle(tree.outputs{leaf})/2,'r.');
% 		end
% 		plot(angle(tree.class(leaves))/2);
% 	axis([ 0,length(leaves)+1,pi/2*[-1,1] ]);
% 	exportfig('U:\matlab\figs\mammography\scatter');

% 	varhst = varhst + hist(tree.var(tree.var~=0),60);

	% compute level of the tree for each node
	tree.level = nan(size(tree.node));
	tree.level(1) = 1;
	for i = 2:length(tree.level)
		tree.level(i) = tree.level(tree.parent(i))+1;
	end
	for i = 1:length(tree.level)
		if tree.var(i)>0
			% branch node - count
			varhst(tree.level(i),tree.var(i)) = ...
				varhst(tree.level(i),tree.var(i)) + 1;
		else
			% leaf node - compute dispersal
			
			% trace all the way back to root, adding dispersal value to every
			% ancestor of node i
			par = tree.parent(i);
			while par>0
				vardisp(tree.level(par),tree.var(par)) = ...
					vardisp(tree.level(par),tree.var(par)) + 1-tree.nodeerr(i);
				vardispcount(tree.level(par),tree.var(par)) = ...
					vardispcount(tree.level(par),tree.var(par)) + 1;
				par = tree.parent(par);
			end
		end
	end
	max_level = max(max_level,max(tree.level)-1);
	
	% fill in co-occurence matrix for neighbouring nodes in tree
	for node = 2:length(tree.node)
		v1 = tree.var(node); % cut variable 1
		if (v1==0), continue; end
		v2 = tree.var(tree.parent(node)); % cut variable 2
		parmat(v1,v2) = parmat(v1,v2) + 1;
	end
end

imagesc(parmat+parmat'); colormap(gray(256)); axis('image');
exportfig([figpath,'neighbours']);

return

varhst = varhst(1:max_level,:);

% normalize dispersal values
vardisp = vardisp./vardispcount;
vardisp = vardisp(1:max_level,:);

figure(2); clf; colormap(gray);
	sbsz = [3,1];
	im1 = varhst./(max(varhst,[],2)*ones(1,D));
	im2 = varhst./(ones(size(varhst,1),1)*max(varhst,[],1));
% 	im1 = vardisp./(max(vardisp,[],2)*ones(1,D));
% 	im2 = vardisp./(ones(max_level,1)*max(vardisp,[],1));
	mysubplot(sbsz,1); imagesc(im1); xlim([0,D+1]);
	mysubplot(sbsz,2); imagesc(im2); xlim([0,D+1]);
	mysubplot(sbsz,3); hold on;
		bar(1:D,sum(varhst,1));
% 		plot(1:D,sum(varhst,1));
		plot([1;1]*(0:6:D)+0.5,[0;1e6]*ones(1,length(0:6:D)),'k:');
		axis([0,D+1,0,max(sum(varhst,1))]);
		

figure; set(gcf,'paperposition',[0,0,20,10]);
	colormap(gray(256));
	clf; 
		imagesc(im1); xlim([0,D+1]); 
		xlabel('Feature'); ylabel('Tree level');
		exportfig([figpath,'histim1']);
	clf; 
		imagesc(im2); xlim([0,D+1]); 
		xlabel('Feature'); ylabel('Tree level');
		exportfig([figpath,'histim2']);
	clf; hold on;
		colormap('default');
		bar(1:D,sum(varhst,1));
% 		plot(1:D,sum(varhst,1));
		plot([1;1]*(0:6:D)+0.5,[0;1e6]*ones(1,length(0:6:D)),'k:');
		axis([0,D+1,0,max(sum(varhst,1))]);
		xlabel('Feature'); ylabel('Count');
		exportfig([figpath,'feathist']);
	clf; hold off;
		h3 = reshape(sum(varhst,1),[6,10]);
		bar3(h3);
		xlabel('Feature'); ylabel('Orientation subband');
		exportfig([figpath,'bar3d']);
	close('all');

return

sbsz = [2,2]; nh = 24;
figure(3); clf;
	mysubplot(sbsz); hist(contrast,nh); title('contrast');
	mysubplot(sbsz); hist(width,nh); title('width');
	mysubplot(sbsz); hist(orientation,nh); title('ori');
	mysubplot(sbsz); hist(squash,nh); title('squash');

