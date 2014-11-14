function classify_retinas(varargin)
% Script for predicting orientation in retinal images and producing stuff
% for the BMVC paper.
%
% Note that for these experiments I'm using regressors trained with 20k
% data points, usually on a 1x1 grid though some will use 3x3

args = u_packargs(varargin,... % the user's input
    '0', ... % non-strict mode
    'job_id',				unixenv('JOB_ID',['pc',datestr(now,'yyyymmddTHHMMSS')]), ...
    'task_id',				unixenv('SGE_TASK_ID',1), ...
	'prediction_type',		unixenv('PREDICTION_TYPE','linear_regression'), ...
	'predictor_job',		unixenv('PREDICTOR_JOB','latest'), ...
	'force_reclassify',		unixenv('FORCE_RECLASSIFY',true), ...
    'dummy',				[]); % keep this on last line for convenience

% if $CUSTOM_ID is defined and not empty then use that instead of $JOB_ID
custom_id = unixenv('CUSTOM_ID');
if ~isempty(custom_id)
	args.job_id = custom_id;
end
job_id = args.job_id;

% whether we need to reclassify or if we can use previous results
f_reclassify = args.force_reclassify;

%% define predictor type and load
pred_type = args.prediction_type;
% pred_type = 'detection'; % random forest classifier
% pred_type = 'logistic_classification'; % logistic classifier
% pred_type = 'orientation'; % random forest regressor
% pred_type = 'linear_regression'; % linear regressor
% pred_type = 'boosted_regression'; % boosted regressor

% find regressor and load it
pred_job	= args.predictor_job;
pred_root	= [asymmetryroot,'data/models/vessel/',pred_type,'/'];
if strcmp(pred_job,'latest')
	% if no specific regressor defined then take the most recent in
	% pred_root
	% Note: this could be risky on Hydra where another node might create a
	% new regressor while these tasks are waiting
	pred_dir	= dir([pred_root,'*']);
	pred_job	= [pred_dir(end).name,'/'];
end
pred_dir	= dir([pred_root,pred_job,'/random_forest*.mat']);

% full filename (including path) to predictor
pred_fname	= [pred_root,pred_job,pred_dir(end).name];
predictor	= u_load(pred_fname);

% get sampling arguments
sampling_args = u_load([pred_root,pred_job,'/sampling_args.mat']);


%% process test images

% retinal root folder
retpath = 'retinograms\drive\test';
retroot = [asymmetryroot,'data/',retpath,'/'];

% path to input images
imgpath = [retpath,'/images_extended'];
imgroot = [asymmetryroot,'data/',imgpath,'/'];

% path where we will save results
outpath = strrep([imgroot,'results/',pred_job,'/'],'\','/');
outpath = strrep(outpath,'//','/');
if ~exist(outpath,'dir')
    f_reclassify = true;
	mkdir(outpath);
	if ~ispc, fileattrib(outpath,'+w','g'); end
end
% path to masks
mskdir = dir([retroot,'mask/*.png']);
lindir = dir([retroot,'1st_manual/*.png']);

% copy parameter list to results dir
if exist([pred_root,pred_job,'args.txt'],'file')
	copyfile([pred_root,pred_job,'args.txt'],...
			 [outpath,'args_',pred_job(1:end-1),'.txt']);
end

% classify set of test images
tic;
if f_reclassify
    classify_image_set(pred_job,imgpath,...
        'forest_dir', sprintf('models/vessel/%s',pred_type), ...
        'mask_dir',	'retinograms\drive\test\vessel_masks\', ...
        'use_nag', ~ispc );
end
T = toc;

%% generate PNG images from saved predictions

% define jazzy colormaps
ori_cmap = [0 0 0; hsv(255)];
mag_cmap = [0 0 0; jet(255)];

% predicted/true orientation over all images (if needed)
ori_est = []; ori_gt = [];

% create folder for error data to go in
if ~exist([outpath,'errors'],'dir')
	mkdir([outpath,'errors']);
end

resdir = dir([outpath,'*_class.mat']);
time_per_image = T/length(resdir);
for i = 1:length(resdir)
	% load predicted output back in
	imgout = load_uint8([outpath,resdir(i).name]);

	% get manual markup
	imgmask = logical(imread([retroot,'mask/',mskdir(i).name]));
	linemask = logical(imread([retroot,'1st_manual/',lindir(i).name]));
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
			fname = [outpath,sprintf('errors/%03d_detection.png',i)];
				imwrite(uint8(1+255*imgout),fname);

		case {'orientation','linear_regression','logistic_regression','boosted_regression'},
			% get angle of output and limit to range [0..pi)
			imgout = mod(angle(imgout)/2,pi);
			imgout_masked = mod(angle(imgout_masked)/2,pi);
			
			% compute error wrt labelled vessels
			load([retroot,'orientations\',sprintf('%02d_ori1.mat',i)]);
			linemask = linemask & ~isnan(gt_ori);
			gt_ori(~linemask) = NaN;
			err_i	= ori_error(imgout(linemask),angle(gt_ori(linemask))/2);
			ori_est	= [ori_est; imgout(linemask)];
			ori_gt	= [ori_gt; angle(gt_ori(linemask))/2];

			% save images to disk
		% 	fname = [outpath,sprintf('%03d_orientation_full.png',i)];
		%		imwrite(uint8(1+255*imgout/pi),cmap,fname);
			fname = [outpath,sprintf('errors/%03d_orientation_masked.png',i)];
				imwrite(uint8(1+255*imgout_masked/pi),ori_cmap,fname);
			fname = [outpath,sprintf('errors/%03d_abs_error.png',i)];
				err_mag = nan(size(linemask)); err_mag(linemask) = abs(err_i);
				imwrite(uint8(1 + 255*err_mag/(pi/2)),mag_cmap,fname);
	end
end

% any final computations
switch pred_type
	case {'detection','logistic_classification'},
		
	case {'orientation','linear_regression','boosted_regression'},
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
						strrep(sampling_args.prediction_type,'_','\_'),...
						sampling_args.rgb_channel,sampling_args.decomp_type,...
						sampling_args.win_size,sampling_args.win_size,...
						es.abs_mean,es.abs_median),...
				'location','southeast');
			title(pred_job);
		graph(1); exportfig([outpath,'errors/cumfreq']);

		save([outpath,'errors/errors.mat'],'err','es','time_per_image');
end

filename = [outpath,'errors/error_stats.txt'];
fid = fopen(filename,'w');
	fprintf(fid,'%s\n',evalc('es'));
	fprintf(fid,'%s\n',evalc('time_per_image'));
fclose(fid);

% create a file with the essential properties of the experiment in its name
filename = sprintf('errors/#%s-%s-%s-%ix%i.txt',...
						sampling_args.prediction_type,...
						sampling_args.rgb_channel,sampling_args.decomp_type,...
						sampling_args.win_size,sampling_args.win_size);
copyfile([outpath,'errors/error_stats.txt'],[outpath,filename]);

