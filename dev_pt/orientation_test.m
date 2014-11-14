%Linop_test
% Script to test 2D line operator.
%Initially set up to threshold the chromosome image, chromos1.tif.
%Make sure it is in the path.
% Edit the source to change the threshold or use other images 

warning('off',...
	      'Images:imshow:magnificationMustBeFitForDockedFigure');

profile off;
% profile clear; profile on;

clear all; clc;


% close all; % get rid of figures from previous runs
iptsetpref('ImshowBorder','tight'); % display preferences

method = 'dt_rf';
% method = 'linop';
% method = 'monogenic';
% method = 'g2d';

src_type = 'grain';
% src_type = 'synthlines'; cropscl = 4;

switch lower(method)
	case 'dt_rf',
		forest_root		= [asymmetryroot,'data/line_orientation_rfs/'];
		forest_dir		= dir([forest_root,'pc*']);
		forest_job		= [forest_dir(end).name,'/'];
%  		forest_job		= 'pc20110131T110913/';
% 		forest_job		= '191934/';
		forest_job		= '242728/';

		forest_folder	= [forest_root,forest_job]
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

imnum = 5;
% theta_err = 0:30:180; theta_err = theta_err(1:end-1); n_tests = 5;
theta_err = 180*rand; n_tests = 1;

theta_err(2:n_tests+1,:) = 0;

% while 1
for i = 1:size(theta_err,2)
	for t = 1:n_tests
		
		switch lower(src_type)
			case 'grain',
				[src,theta] = create_grain(15,theta_err(1,i));
				src = double(src);
				ori_gt = theta*ones(size(src));
				ori_map = logical(ones(size(src)));

			case 'synthlines',
				bgpath = 'U:\projects\mammography\data\synthetic_lines2\';
				bgdir = dir([bgpath,'*.mat']);
				load([bgpath,bgdir(imnum).name]);

				lblpath = 'U:\projects\mammography\data\synthetic_lines2\labels\';
				lbldir = dir([lblpath,'*.mat']);
				load([lblpath,lbldir(imnum).name]); 

				imsz = floor(size(test_image)/cropscl);
				src = test_image(1:imsz(2),1:imsz(1));
				ori_gt = label_orientation(1:imsz(2),1:imsz(1));
				ori_map = (label(1:imsz(2),1:imsz(1))==1);
				imnum = imnum+1;

			otherwise,
			% src = imread('102_2.tif'); % load the image; 
		end

		tic;
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
				figure(10); hist(line_scale(:));

			case 'linop',
				num_angles = 24; % number of orientations for detecting lines
				min_scale = 5; % smallest scale for detecting linear segments;
				max_scale = 15; % largest scale for detecting linear segments;

				[response, orientations] = line_operator_conv(src, num_angles,min_scale, max_scale, 'degrees');
				line_im = non_maximal_supp(response, orientations);

			otherwise,
				% do something
		end
		toc;

		ori_img = orientations;
		ori_img(~ori_map) = NaN;
		ori_err = mean(abs(mb_mod(orientations(ori_map)-ori_gt(ori_map),180)))
		
		theta_err(1+t,i) = ori_err;
		if size(theta_err,1)>2 || size(theta_err,2)>1
			continue;
		end

		figure(1);
			spsz = [2,2];
			mysubplot(spsz,1); imagesc(src); axis('image','ij','off'); title('input'); % display the image
			mysubplot(spsz,2); imagesc(-response); axis('image','ij','off'); title('line response');% display the response image
			mysubplot(spsz,3); imshow(ori_img/180); title('orientation');
			mysubplot(spsz,4); imshow(line_im); title('line\_im');

		figure(2);
			subplot(2,2,1); show_ori_image(ori_img);
			subplot(2,2,3); show_ori_image(ori_gt);
			subplot(2,2,2); show_ori_image;
			subplot(2,2,4); show_ori_image;

		if strcmp(src_type,'grain')
			figure(3); clf; hold on;
				h = histc(orientations(:),0:2:180)';
				h = [h(1:end-1) h(1:end-1)];
				t = linspace(0,360,length(h)+1); t = t(1:end-1); t = t+(t(2)-t(1))/2;
				plot(	theta([1,1]),10e6*[0,1],'r:',...
							theta([1,1])+180,10e6*[0,1],'r:');
				plot(	mod(theta([1,1])+90,360),10e6*[0,1],':',...
							mod(theta([1,1])+270,360),10e6*[0,1],':',...
							'color',0.7*[1,1,1]);
				plot(	t,h,'b-');
				axis([0,360,0,max(h)*1.05]);
		end
	end % tests	
% 	break
end	% theta

profstat = profile('status');
if strcmp(profstat.ProfilerStatus,'on'), profile report; end

if size(theta_err,2)>1
	figure; plot(theta_err(1,:),theta_err(2:n_tests+1,:));
end
