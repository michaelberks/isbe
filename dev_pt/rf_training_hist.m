function [in_hist,out_hist,out_hist_w] = rf_training_hist(rf)
% generate histograms of training data and RF output
% histogram bins correspond to angles 0..179 degrees

in_hist = zeros(180,1);
out_hist = zeros(180,1);
out_hist_w = zeros(180,1);
for t = 1:length(rf.trees)
	filename = [rf.tree_root,rf.tree_dir,rf.trees{t}];
	load(filename);
	
	% update input histogram
	if isfield(tree,'y_hist')
		in_hist = in_hist + tree.y_hist;
	end
	
	% use only leaf nodes for output histogram
	leaves = find(all(tree.children==0,2));

	% convert from complex value to degrees
	if ~isreal(tree.class)
		tree.class = angle(tree.class)/2*(180/pi)+90;
	end
	
	% update output histogram
	hst = histc(tree.class(leaves),-0.5:180.5);
	hst(end) = []; % last bin corresponds to theta=180.5 exactly
	hst(1) = hst(1)+hst(end); % copy entries from theta=180 to theta=0 	
	hst(end) = []; % delete theta=180 bin
	out_hist = out_hist + hst;

	% update weighted output histogram if data available
	if isfield(tree,'nodesize')
		hst(:) = 0; % clear unweighted histogram
		for ileaf = 1:length(leaves)
			leaf = leaves(ileaf);
			bin = round(tree.class(leaf));
			if bin==180, bin = 0; end % bins go from 0..179
			hst(bin+1) = hst(bin+1) + tree.nodesize(leaf); % theta=0 => bin 1
		end
	end
	out_hist_w = out_hist_w + hst;
end

% normalize histograms
if all(in_hist==0),	in_hist(:) = NaN;
else,								in_hist = in_hist/sum(in_hist);
end
out_hist = out_hist/sum(out_hist);
out_hist_w = out_hist_w/sum(out_hist_w);

if nargout==0,
	% display results
	figure;
		subplot(3,1,1); plot(0:179,in_hist); title('Input dist.');
		subplot(3,1,2); plot(0:179,out_hist); title('Output dist.');
		subplot(3,1,3); plot(0:179,out_hist_w); title('Output dist. (weighted)');
	clear in_hist out_hist out_hist_w;
end