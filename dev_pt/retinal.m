% define image path and get list of PNG images
imgpath	= 'U:\data\retinal\drive\test\';
imgdir	= dir([imgpath,'images\*.png']);
mskdir	= dir([imgpath,'mask\*.png']);
lindir	= dir([imgpath,'1st_manual\*.png']);

% load most recently created forest if none specified
clear tree;
forest_root		= [asymmetryroot,'data/line_orientation_rfs/'];
forest_dir		= dir([forest_root,'pc*']);
forest_job		= [forest_dir(end).name,'/'];
forest_job		= ['pc20110307T145627/']; % 40x16000
forest_dir		= dir([forest_root,forest_job,'/random_forest*.mat']);
forest_fname	= [forest_root,forest_job,forest_dir(end).name];
forest			= u_load(forest_fname);

n_trees = length(forest.trees);
n_trees = 1;

% get sampling arguments
sampling_args = u_load([forest_root,forest_job,'/sampling_args.mat']);

%% generate training data
fprintf('Generating training data...');
sampling_args = struct(...
	'saved_data_dir','A:\data\synthetic_data\drive\dt\training\',...
    'num_samples',16000,...
    'task_id',[],...
	'pts_per_image',4000, ...
    'feature_type', 'conj',...
    'win_size', 1);
[X,y] = sample_saved_dt_line_data(sampling_args);
fprintf('done\n');

%% build regressors
fprintf('Building...');
for t = 1:n_trees
% 	load([forest.tree_root,forest.tree_dir,sprintf('traindata%04d.mat',t)]);

	% linear regressor
	fprintf('linear...');
	linreg = pt_lin_reg_train(X,y);
end
fprintf('done\n');

% sampling arguments for DTCWT coefficients
sampling_args = get_substructure(...
	sampling_args,{'win_size','feature_type'});

for i = 1:1
	% load image and take green channel
	img = imread([imgpath,'images/',imgdir(i).name]);
	img = img(:,:,2);
	
	% get manual markup
	imgmask = logical(imread([imgpath,'mask/',mskdir(i).name]));
	linemask = logical(imread([imgpath,'1st_manual/',lindir(i).name]));
	
	% display
	figure(1); clf; colormap(gray(256));
		image(img); axis('image');

	% matrices of row and column coefficients
	[cc,rr] = meshgrid(1:size(img,2),1:size(img,1));
	
	% compute DTCWT coefficients
	dt = compute_dual_tree(img,5,0);
	
	% output image
	imgout = nan(size(img));

	% compute coefficients in blocks
	blocksz = 128;
	nblocks = ceil(size(img)/blocksz);
	
	rows = 1:blocksz;
	rows(rows>size(img,1)) = [];
	for bx = 1:nblocks(2)

		cols = 1:blocksz;
		cols(cols>size(img,2)) = [];
		for by = 1:nblocks(1)
			% sampled pixels
			inds = sub2ind(size(img),rr(rows,cols),cc(rows,cols));
			[rinds,cinds] = ind2sub(size(img),inds(:));

			% get these coefficients
			Xtest = sample_dt_data(dt,rinds,cinds,sampling_args);
	
			% predict orientations from coefficients
			y_linear = linear_regressor_predict(linreg,Xtest);
	
			% copy to output image
			imgout(inds) = angle(y_linear)/2;
			
			% next column block
			cols = cols+blocksz;
			cols(cols>size(img,2)) = [];
		end
		
		% next row block
		rows = rows+blocksz;
		rows(rows>size(img,1)) = [];
	end
	
	% mask out background pixels
	imgout(~imgmask) = NaN;
	imgout2 = imgout;
	imgout2(~linemask) = NaN;
	
	% display
	cmap = hsv(256); cmap(1,:) = 0;
	figure(2); clf; colormap(cmap);
	subplot(1,2,1);
		imagesc([imgout; linspace(-pi,pi,size(imgout,2))/2]); axis('image');
	subplot(1,2,2);
		imagesc([imgout2; linspace(-pi,pi,size(imgout2,2))/2]); axis('image');
end

