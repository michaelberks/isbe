warning('off',...
	      'Images:imshow:magnificationMustBeFitForDockedFigure');

profile off;
% profile clear; profile on;

%% initialization
clear all; clc;

% seed random numbers with fixed value
rand('twister',5489);

% close all; % get rid of figures from previous runs
iptsetpref('ImshowBorder','tight'); % display preferences

% method = 'dt_rf';
% method = 'linop';
% method = 'monogenic';
method = 'g2d';

%% load model
switch lower(method)
	case 'dt_rf',
		% use local models...
		forest_root		= [asymmetryroot,'data/line_orientation_rfs/'];
		forest_dir		= dir([forest_root,'pc*']);
		forest_job		= [forest_dir(end).name,'/'];

		% ...or use one of Mike's
		forest_job		= '191934/';

		forest_folder	= [forest_root,forest_job];
		forest_fname	= [forest_folder,'/random_forest.mat'];
		forest				= u_load(forest_fname);
		
		sampling_args = u_load([forest_folder,'/sampling_args.mat']);
		sampling_args_c = get_substructure( sampling_args,...
			{'num_levels','feature_shape','feature_type','do_max','rotate',...
			 'win_size','use_nag'} );
		
		forest.tree_root = [asymmetryroot,'data/line_orientation_rfs/'];
		
	otherwise,
		% do something
end

%% load a selection of backgrounds
n_bgs = 20;
bgroot = [asymmetryroot,'\data\synthetic_backgrounds\real512\test\'];
bgdir = dir([bgroot,'bg*.mat']);
for ibg = 1:n_bgs
	fname = bgdir(ceil(rand*length(bgdir))).name;
	tmp = load([bgroot,fname]);
	if isfield(tmp,'bg')
		bg{1,ibg} = tmp.bg(1:128,1:128);
	elseif isfield(tmp,'sampled_window')
		bg{1,ibg} = tmp.sampled_window(1:128,1:128);
	end
	
	% and a rotated version, too
	bg{2,ibg} = bg{1,ibg}';
end

%% generate synthetic line images
theta = 0:30:180;

% number of images per theta
n_per_theta = 10;

img = {}; lbl = {};
for itheta = 1:length(theta)
	for iimg = 1:n_per_theta
		% choose background at random
		ibg = ceil(rand*2*n_bgs); % two per background
		img{itheta,iimg} = bg{ibg};

		% generate line image and its label
		[fg,label] = create_ellipse_bar(5, 4, theta(itheta), size(img{itheta},1),size(img{itheta},2), 64,64);	

		% superimpose a line
		img{itheta,iimg} = img{itheta,iimg} + fg;
		lbl{itheta,iimg} = label;
		ori{itheta,iimg} = theta(itheta)*label;
	end
end

% return

%% process the images
errmat = zeros(n_per_theta,length(theta));
n_remaining = length(theta)*n_per_theta;
n_done = 0;

% create waitbar
if exist('wb','var'), close(wb); end
msg = sprintf('Estimated time remaining = ???');
wb = waitbar(0,msg,'position',[975 52 270 56.25]);
tic;
for itheta = 1:length(theta)
	for iimg = 1:n_per_theta
	
		src = img{itheta,iimg};
		ori_map = lbl{itheta,iimg};
		ori_gt = ori{itheta,iimg};
		
		switch lower(method)
			case 'dt_rf',
				[orientations] = classify_image(...
						'image_in', src, ...
						'forest', forest,...
						'sampling_args', sampling_args_c,...
						'forest_type', 'regression',...
						'decomp_type', 'dt');

				response = ori_map;
				line_im = zeros(size(src));

			case 'monogenic',
				[image_out,phase_cong,orientations,local_ori] = ...
					monogenic_phase_cong(src,3,4,2,0.65);

				response = zeros(size(src));
				line_im = zeros(size(src));

			case 'g2d',
				[orientations, line_im, line_scale, response] = karssemeijer_line_detection(double(src));
				orientations = mod(orientations+360,180);
% 				figure(10); hist(line_scale(:));

			case 'linop',
				num_angles = 24; % number of orientations for detecting lines
				min_scale = 5; % smallest scale for detecting linear segments;
				max_scale = 15; % largest scale for detecting linear segments;

				[response, orientations] = line_operator_conv(src, num_angles,min_scale, max_scale, 'degrees');
				line_im = non_maximal_supp(response, orientations);

			otherwise,
				% do something
		end

		ori_err = mean(abs(mb_mod(orientations(ori_map)-ori_gt(ori_map),180)));
		
		errmat(iimg,itheta) = ori_err;

		% update waitbar
		n_done = n_done+1; n_remaining = n_remaining-1;
		t = toc * n_remaining/n_done;
		mm = floor(t/60); ss = t - mm*60;
		msg = sprintf('Estimated time remaining = %.0fm %.0fs',mm,ss);
		waitbar(n_done/(n_done+n_remaining),wb,msg);
	end % for iimg
end	% for itheta
close(wb);

profstat = profile('status');
if strcmp(profstat.ProfilerStatus,'on'), profile report; end

figure(1); clf; hold on;
	errorbar(theta,mean(errmat),std(errmat));
	axis([-5,185,0,90]);
