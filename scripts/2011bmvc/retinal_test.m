% Script for predicting orientation in retinal images and producing stuff
% for the BMVC paper.
%
% Note that for these experiments I'm using regressors trained with 20k
% data points, usually on a 1x1 grid though some will use 3x3

clear variables;

%% define predictor type and load
% pred_type = 'detection'; % random forest classifier
% pred_type = 'logistic_classification'; % logistic classifier
pred_type = 'orientation'; % random forest regressor
% pred_type = 'linear_regression'; % linear regressor
% pred_type = 'boosted_regression'; % boosted regressor
% pred_type = 'analytic';

asymmetryroot = 'Z:/asym/';

switch pred_type,
	case 'orientation',
% 		pred_job	= 'pc20110307T145627/'; % 40x16000 forest
% 		pred_job	= '313710/'; % 40x16000 forest
% 		pred_job	= 'dt';
% 		pred_job	= '313846/'; % 200k mono, 3x3 windows
%       pred_job    = '48947'; % clover (with 5sigma support)
%       pred_job    = '48949'; %200k momo, 1x1 windows
        pred_job    = '48948'; %200k g1d, 3x3 windows
		
	case 'linear_regression',
% 		pred_job	= '313814/'; % 200k dt, generated
%  		pred_job	= '313804/'; % 200k haar ('right' sigmas)
% 		pred_job	= '313803/'; % 200k clover ('right' sigmas)
% 		pred_job	= '313836/'; % 200k mono

	case 'boosted_regression',
% 		pred_job	= '313813/'; % 200k dt, generated
%  		pred_job	= '313820/'; % 200k haar ('right' sigmas)
% 		pred_job	= '313829/'; % 200k clover ('right' sigmas)
% 		pred_job	= '313838/'; % 200k mono

	case 'analytic',
% 		pred_job	= 'clover';
% 		pred_job	= 'mono';
        pred_job    = 'g1d';
end

% find regressor
pred_root	= [asymmetryroot,'data/models/vessel/',pred_type,'/'];
% if no specific regressor defined then take the most recent
if ~exist('pred_job','var'), 
	pred_dir = dir([pred_root,'pc*']);
	pred_job = [pred_dir(end).name,'/']; 
end

% define and create output folder
retpath = 'retinograms\drive\test';
retroot = [asymmetryroot,'data/',retpath,'/'];
imgpath = [retpath,'/images_extended'];
imgroot = [asymmetryroot,'data/',imgpath,'/'];
outpath = strrep([imgroot,'results/',pred_job,'/'],'\','/');
outpath = strrep(outpath,'//','/');
if ~exist(outpath,'dir'), mkdir(outpath); end

mskdir = dir([retroot,'foveal_masks/*.mat']);
lindir = dir([retroot,'vessel_masks/*.mat']);

sampling_args = struct(...
	'decomp_type','',...
	'rgb_channel','rgb',...
	'win_size',3,...
	'prediction_type',pred_type);

T = -1;

% if predictor folder exists then load it and its sampling_args
% if exist([pred_root,pred_job],'dir')
% 	pred_dir	= dir([pred_root,pred_job,'/random_forest*.mat']);
% 	pred_fname	= [pred_root,pred_job,pred_dir(end).name];
% 	predictor	= u_load(pred_fname);
% 
% 	% get sampling arguments
% 	sampling_args = u_load([pred_root,pred_job,'/sampling_args.mat']);
% 
% %% process test images
% 
% 	% copy parameter list to results dir
% 	copyfile([pred_root,pred_job,'args.txt'],outpath);
% 
% 	% classify set of test images
% 	resdir = dir([outpath,'*_class.mat']);
% 	if ~strcmp(pred_type,'orientation')
% 		tic;
% 		classify_image_set(pred_job,imgpath,...
% 			'forest_dir', sprintf('models/vessel/%s',pred_type), ...
% 			'mask_dir',	'retinograms\drive\test\vessel_masks\', ...
% 			'use_nag', ~ispc );
% 		T = toc;
% 	end
% end

%% generate PNG images from saved predictions

% define jazzy colormaps
ori_cmap = [0 0 0; hsv(255)];
mag_cmap = [0 0 0; jet(255)];

% predicted/true orientation over all images (if needed)
ori_est = []; ori_gt = [];

if ~isfield(sampling_args,'rgb_channel')
	sampling_args.rgb_channel = 'rgb';
end

resdir = dir([outpath,'*_class.mat']);
if isempty(resdir)
	resdir = dir([outpath,'*_ori.mat']);
end
time_per_image = T/length(resdir);
for i = 1:length(resdir)
	% load predicted output back in
	vars = whos('-file',[outpath,resdir(i).name]);
	var_names = {vars(:).name};
	if any(strcmp(var_names,'scaling'))
		imgout = load_uint8([outpath,resdir(i).name]);
	elseif any(strcmp(var_names,'ori_map'))
		load([outpath,resdir(i).name]);
		imgout = ori_map;
	end

	% get manual markup
	imgmask = u_load([retroot,'foveal_masks/',mskdir(i).name]);
	linemask = u_load([retroot,'vessel_masks/',lindir(i).name]);
	linemask(~imgmask) = false;

	% mask out background pixels and non-line pixels
	imgout(~imgmask) = NaN;
	imgout_masked = imgout;
	imgout_masked(~linemask) = NaN;

	% process raw output depending on regression/classification task
	switch pred_type
		case {'detection','logistic_classification'},
			% convert to three colour format
			% y = TP; k = TN; g = FN; r = FP;
			imgout(:,:,2) = linemask;
			imgout(:,:,3) = 0;

			% save image to disk
			fname = [outpath,sprintf('%03d_detection.png',i)];
				imwrite(uint8(1+255*imgout),fname);

		case {'analytic','orientation','linear_regression',...
				'logistic_regression','boosted_regression'},
			if strcmp(pred_type,'analytic')
				% angle needs comverting to complex
				imgout = exp(complex(0,imgout*2));
				imgout_masked = exp(complex(0,imgout_masked*2));
			elseif strcmp(pred_type,'orientation') && length(pred_job) > 5
				% angle needs doubling for 'orientation'
				imgout = exp(complex(0,angle(imgout)*2));
				imgout_masked = exp(complex(0,angle(imgout_masked)*2));
			end				

			% compute error wrt labelled vessels
			load([retroot,'orientations\',sprintf('%02d_ori1.mat',i)]);
			linemask = linemask & ~isnan(gt_ori);
			gt_ori(~linemask) = NaN;
			err_i	= ori_error(imgout(linemask),gt_ori(linemask));
			ori_est	= [ori_est; imgout(linemask)];
			ori_gt	= [ori_gt; gt_ori(linemask)];

			% save images to disk
		% 	fname = [outpath,sprintf('%03d_orientation_full.png',i)];
		%		imwrite(uint8(1+255*imgout/pi),cmap,fname);
			fname = [outpath,sprintf('%03d_orientation_masked.png',i)];
				imwrite(uint8(1+255*mod(angle(imgout_masked)/2,pi)/pi),ori_cmap,fname);
			fname = [outpath,sprintf('%03d_abs_error.png',i)];
				err_mag = nan(size(linemask)); err_mag(linemask) = abs(err_i);
				imwrite(uint8(1 + 255*err_mag/(pi/2)),mag_cmap,fname);
	end
end

% any final computations
switch pred_type
	case {'detection','logistic_classification'},
		
	case {'analytic','orientation','linear_regression','boosted_regression'},
		% compute errors, draw graph and save
		[err,es] = ori_error(ori_est,ori_gt);
		figure(1); clf; hold on;
			h = []; % handles to plots
			plot([0,90],0.5*[1,1],':','color',0.8*[1,1,1]); % median
			plot([0,90],0.9*[1,1],':','color',0.8*[1,1,1]); % 90%ile
			h = plot(es.abs_percentiles,0.01:0.01:1);
			axis([0,90,0,1]);
			xlabel('Error (degrees)'); ylabel('Cumulative Frequency');
			legend(h,sprintf('%s (%s, %s, %ix%i): mean %0.2f, med %0.2f',...
						strrep(pred_type,'_','\_'),...
						sampling_args.rgb_channel,sampling_args.decomp_type,...
						sampling_args.win_size,sampling_args.win_size,...
						es.abs_mean,es.abs_median),...
				'location','southeast');
			title(pred_job);
		graph(1); exportfig([outpath,'cumfreq']);

		save([outpath,'errors.mat'],'err','es','time_per_image');
end

filename = [outpath,'error_stats.txt'];
fid = fopen(filename,'w');
	fprintf(fid,'%s\n',evalc('es'));
	fprintf(fid,'%s\n',evalc('time_per_image'));
fclose(fid);

% create a file with the essential properties of the experiment in its name
filename = sprintf('#%s-%s-%s-%ix%i.txt',...
						pred_type,...
						sampling_args.rgb_channel,sampling_args.decomp_type,...
						sampling_args.win_size,sampling_args.win_size);
copyfile([outpath,'error_stats.txt'],[outpath,filename]);

