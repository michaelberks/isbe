function tree = show_tree_output(tree,max_levels,min_radius)
%SHOW_TREE_OUTPUT Displays the outputs of a tree with orientation output as
% a 'wheel' consisting of concentric circles (one per level), each showing
% the output of all leave nodes as a coloured marker.
%   [] = show_tree_output(tree,max_levels,min_radius)
%
% Inputs:
%      tree - Tree structure that predicts orientation output
%
%      max_levels - Greatest number of levels to show
%
%      min_radius - Inner radius of the 'wheel' (0 < min_radius < 1,
%      default 0.2)
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 17-Feb-2011
% Author: Phil Tresadern 
% Email : philip.tresadern@manchester.ac.uk 
% Phone : +44 (0)161 275 5114 
% Copyright: (C) University of Manchester 

if nargin<1 || isempty(tree)
	% load most recently created tree if none specified
	forest_root		= [asymmetryroot,'data/line_orientation_rfs/'];
	forest_dir		= dir([forest_root,'pc*']);
	forest_job		= [forest_dir(end).name,'/'];
	forest_dir		= dir([forest_root,forest_job,'/random_forest*.mat']);
	forest_fname	= [forest_root,forest_job,forest_dir(end).name];
	forest				= u_load(forest_fname);
	tree					= u_load([forest.tree_root,forest.tree_dir,forest.trees{1}]);
end

% set default max_levels
if nargin<2 || max_levels<1
	max_levels = inf;
end

% use default min_radius
if nargin<3 || (min_radius<0) || (min_radius>1)
	min_radius = 0.2; 
end


% compute level of each node
treelevel = nan(size(tree.node));
treelevel(1) = 0;
for i = 2:length(treelevel)
	treelevel(i) = treelevel(tree.parent(i))+1;
end

% get total number of levels up to some maximum (default=inf)
n_levels = min(max(treelevel),max_levels);

% get radii of concentric circles
radii = linspace(min_radius,1,n_levels+1);

% set up the figure
set(gcf,'windowstyle','normal'); 
clf; hold on;

% draw the circle boundaries
t = linspace(0,2*pi,101);
for i = 1:length(radii)
	plot(radii(i)*cos(t),radii(i)*sin(t),':','color',0.8*[1,1,1]);
end

cols = {'r-','b-'};

% find leaf nodes (ish)
leaves = find(tree.var==0);

% plot each leaf at every level
for L = 1:length(leaves)
	node = leaves(L);

	% get output value and compute line direction
	theta = angle(tree.class(node));
	Tdir = [cos(theta); sin(theta)];
	while tree.parent(node)>0
		if (treelevel(node)<=n_levels)
			% get the level and compute the (x,y) coords of the line endpoints
			lvl = treelevel(node);
			x1 = radii(lvl)*Tdir; x2 = radii(lvl+1)*Tdir;

			% find whether this is the left or right child of the parent and set
			% the colour accordingly
			ichild = find(tree.children(tree.parent(node),:)==node);
			linestyle = cols{ichild};

			% draw the line for this leaf at this level
			plot([x1(1),x2(1)],[x1(2),x2(2)],linestyle);
		end
		
		% go to parent node and repeat
		node = tree.parent(node);
	end
end

% reshape the axis
axis(1.1*[-1,1,-1,1],'square');





