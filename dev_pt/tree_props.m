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
	
	display(forest_fname);
end

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
% 	plot(hist(tree.level,max(tree.level)));
	
	levelsize = zeros(1,max(tree.level));
	for i = 1:max(tree.level)
		inds = find(tree.level==i);
		levelsize(i) = mean(tree.nodesize(inds));
	end
	
% 	set(figure(1),'windowstyle','docked'); clf; hold on;
% 	ls_ideal = tree.nodesize(1)./(2.^[0:max(tree.level)-1]);
% 	plot(	0:max(tree.level)-1,log(ls_ideal)/log(2),'r-',...
% 				0:max(tree.level)-1,log(levelsize)/log(2),'b:');
% 	ylim([0,log(tree.nodesize(1))/log(2)]);
% 	return
% 	
% 	plot(0:length(levelsize)-1,log(tree.nodesize(1)./levelsize)/log(2));
% 	plot([0,1e6],[0,1e6],'r-');
% 	xmax = max(tree.level); ymax = log(tree.nodesize(1))/log(2);
% 	axis([0,xmax,0,ymax]);
% 	return
	
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

varhst = varhst(1:max_level,:);

% normalize dispersal values
vardisp = vardisp./vardispcount;
vardisp = vardisp(1:max_level,:);

% get marginal histogram
marginal = sum(varhst,1);

figure(2); clf; colormap(gray);
	sbsz = [3,1];
	im1 = varhst./(max(varhst,[],2)*ones(1,D));
	im2 = varhst./(ones(size(varhst,1),1)*max(varhst,[],1));
% 	im1 = vardisp./(max(vardisp,[],2)*ones(1,D));
% 	im2 = vardisp./(ones(max_level,1)*max(vardisp,[],1));
	mysubplot(sbsz,1); imagesc(im1); xlim([0,D+1]);
	mysubplot(sbsz,2); imagesc(im2); xlim([0,D+1]);
	mysubplot(sbsz,3); hold on;
		bar(1:D,marginal);
% 		plot(1:D,marginal);
		plot([1;1]*(0:6:D)+0.5,[0;1e6]*ones(1,length(0:6:D)),'k:');
		axis([0,D+1,0,max(marginal)]);

return
		
figure(1); set(gcf,'paperposition',[0,0,20,10]);
% 	colormap(gray(256));
% 	clf; 
% 		imagesc(im1); xlim([0,D+1]); 
% 		xlabel('Feature'); ylabel('Tree level');
% 		exportfig([figpath,'histim1']);
% 	clf; 
% 		imagesc(im2); xlim([0,D+1]); 
% 		xlabel('Feature'); ylabel('Tree level');
% 		exportfig([figpath,'histim2']);
	clf; hold on;
		colormap('default');
		bar(1:D,marginal);
% 		plot(1:D,marginal);
		plot([1;1]*(0:6:D)+0.5,[0;1e6]*ones(1,length(0:6:D)),'k:');
		axis([0,D+1,0,max(marginal)]);
		xlabel('Feature'); ylabel('Count');
% 		exportfig([figpath,'feathist']);
% 	clf; hold off;
% 		h3 = reshape(marginal,[6,10]);
% 		bar3(h3);
% 		xlabel('Feature'); ylabel('Orientation subband');
% 		exportfig([figpath,'bar3d']);
% 	close('all');

return

[ignore,inds] = sort(marginal) % sort in order of popularity

% look at what options were available for unpopular variables
inds = 1:D;
pairmat = zeros(D,D);
for i = 1:length(inds)
	if marginal(inds(i))==0, continue; end
	nodeopts = tree.nodeopts(find(tree.var==inds(i)),1:2:end);
	pairmat(i,:) = pairmat(i,:) + hist(nodeopts(:),1:D);
	pairmat(i,i) = 0;
	pairmat(i,:) = pairmat(i,:)/sum(pairmat(i,:));
end

figure(2); clf; colormap(gray(256));
	imagesc(pairmat); axis('image');
	ylabel('Selected'); xlabel('Other candidates');
	exportfig([figpath,'pairmat']);

figure(2); clf; colormap(gray(256));
	imagesc(parmat+parmat'); axis('image');
	xlabel('Node'); ylabel('Neighbour');
	exportfig([figpath,'neighbours']);

sbsz = [2,2]; nh = 24;
figure(3); clf;
	mysubplot(sbsz); hist(contrast,nh); title('contrast');
	mysubplot(sbsz); hist(width,nh); title('width');
	mysubplot(sbsz); hist(orientation,nh); title('ori');
	mysubplot(sbsz); hist(squash,nh); title('squash');

