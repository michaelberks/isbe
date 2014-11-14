% Script for predicting orientation in synthetic line images and producing
% stuff for the BMVC paper.
%
% Note that for these experiments I'm using regressors trained with 20k
% data points, usually on a 1x1 grid though some will use 3x3

close all hidden;
clear variables;

%% define predictor type and load
% pred_type = 'detection'; % random forest classifier
% pred_type = 'logistic_classification'; % logistic classifier
% pred_type = 'orientation'; % random forest regressor
% pred_type = 'linear_regression'; % linear regressor
% pred_type = 'boosted_regression'; % boosted regressor
pred_type = 'analytic'; % no regressor

results_ind = '';
switch pred_type,
	case 'orientation',
% 		pred_job	= '286712/'; % dt, optimal vars
% 		pred_job	= '287361/'; % g2d
% 		pred_job	= '313715/'; % haar
% 		pred_job	= '10762/'; % clover
% 		pred_job	= '10766/'; % monogenic

% 		pred_job	= 'pc20110307T145627/'; % 40x16000 forest
% 		pred_job	= '306060/'; results_ind = '005'; % dt, random vars
		
	case 'linear_regression',
% 		pred_job	= '313751/'; % 200k dt, generated
%  		pred_job	= '313759/'; % 200k haar ('right' sigmas)
% 		pred_job	= '313801/'; % 200k clover ('right' sigmas)
% 		pred_job	= '313850/'; % 200k monogenic

		% for speed tests
% 		pred_job	= 'pc20110421T220756/'; % g2d
% 		pred_job	= 'pc20110421T201954/'; % clover
% 		pred_job	= 'pc20110421T202049/'; % haar
% 		pred_job	= 'pc20110421T202816/'; % dt
% 		pred_job	= 'pc20110421T212449/'; % mono
		
% 		pred_job	= '313756/'; % 200k dt, sampled
% 		pred_job	= '313757/'; % 40k dt, generated
% 		pred_job	= '313758/'; % 40k clover ('wrong' sigmas)
% 		pred_job	= 'pc20110420T223839/'; % 40k clover

	case 'boosted_regression',
% 		pred_job	= 'pc20110421T094247/'; % 200k clover
% 		pred_job	= '313763/'; % 200k haar
% 		pred_job	= '313764/'; % 200k dt, sampled
% 		pred_job	= '313765/'; % 200k clover
		pred_job	= '313851/'; % 200k monogenic
		
	case 'analytic',
% 		pred_job	= 'clover';
		pred_job	= 'mono2';
end

% default sampling arguments
sampling_args = struct(...
	'decomp_type','',...
	'rgb_channel','',...
	'win_size',-1,...
	'prediction_type',pred_type);

if ~strcmp(pred_type,'analytic')
	% find regressor and load it
	pred_root	= [asymmetryroot,'data/models/mammograms/',pred_type,'/'];
	% if no specific regressor defined then take the most recent
	if ~exist('pred_job','var'), 
		pred_dir = dir([pred_root,'pc*']);
		pred_job = [pred_dir(end).name,'/']; 
	end

	pred_path	= [pred_root,pred_job];
	if exist(pred_path,'dir')
		% if we have the regressor then load it
		pred_dir	= dir([pred_path,'/random_forest*.mat']);
		pred_fname	= [pred_path,pred_dir(end).name];
		predictor	= u_load(pred_fname);

		% get sampling arguments
		sampling_args = u_load([pred_path,'/sampling_args.mat']);
	end

	% backward compatibility
	if ~isfield(sampling_args,'prediction_type')
		sampling_args.prediction_type = sampling_args.detection_type;
	end
	if ~isfield(sampling_args,'decomp_type')
		sampling_args.decomp_type = 'dt';
	end
end

%% process test images

% define and create output folder
retpath = 'synthetic_lines/real512';
retroot = [asymmetryroot,'data/',retpath,'/'];
imgpath = [retpath];
imgroot = [asymmetryroot,'data/',imgpath,'/'];
outpath = strrep([imgroot,'results/',pred_job,'/',results_ind],'\','/');
outpath = strrep(outpath,'//','/');
if ~exist(outpath,'dir'), mkdir(outpath); end

if exist('pred_path','var') & exist(pred_path,'dir')
	% copy parameter list to results dir
	if exist([pred_path,'args.txt'],'file')
		copyfile([pred_path,'args.txt'],outpath);
	end
end

% get list of masks
mskdir = dir([retroot,'masks/*.mat']);
lbldir = dir([retroot,'labels/*.mat']);
	
% view = '001';
% view = '1';
view = [];

T = -1;
resdir = dir([outpath,'*_class.mat']);
if isempty(resdir)
	resdir = dir([outpath,'*_ori.mat']);
end
if length(resdir)<100 && ~strcmp(pred_type,'analytic')
	% classify set of test images
	tic;
	classify_image_set(pred_job,imgpath,...
		'view',view,...
		'forest_dir', sprintf('models/mammograms/%s',pred_type), ...
		'mask_dir',	'synthetic_lines/real512/masks', ...
		'use_nag', ~ispc);
	T = toc;
end


%% generate PNG images from saved predictions

% define jazzy colormaps
ori_cmap = [0 0 0; hsv(255)];
mag_cmap = [0 0 0; jet(255)];

% predicted/true orientation over all images (if needed)
ori_est = []; ori_gt = [];

resdir = dir([outpath,'*_class.mat']);
if isempty(resdir)
	resdir = dir([outpath,'*_ori.mat']);
end
time_per_image = T/length(resdir);
tb = timebar('limit',length(resdir),'title','Converting images');
for i = 1:length(resdir)
	% load prediced output back in
	vars = whos('-file',[outpath,resdir(i).name]);
	var_names = {vars(:).name};
	if any(strcmp(var_names,'scaling'))
		imgout = load_uint8([outpath,resdir(i).name]);
	elseif any(strcmp(var_names,'orientation_image'))
		load([outpath,resdir(i).name]);
		imgout = orientation_image;
	end

	% get manual markup
	[p,f,e] = fileparts(mskdir(i).name);
	switch e
		case '.png',
			linemask = logical(imread([retroot,'masks/',mskdir(i).name]));
		case '.mat',
			linemask = u_load([retroot,'masks/',mskdir(i).name]);
	end

	% mask out background pixels and non-line pixels
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
				imgout = exp(complex(0,2*imgout/180*pi));
				imgout_masked = exp(complex(0,2*imgout_masked/180*pi));
			elseif strcmp(pred_type,'orientation')
				% angle needs doubling for 'orientation'
				imgout = exp(complex(0,angle(imgout)*2));
				imgout_masked = exp(complex(0,angle(imgout_masked)*2));
			end				
% 			figure(3); image(uint8(1+255*mod(angle(imgout_masked)*0.5,pi)/pi));
			
			% compute error wrt labelled vessels
			load([retroot,'labels\',lbldir(i).name]);
			gt_ori	= complex(	cosd(2*label_orientation),...
								sind(2*label_orientation) );
			linemask = linemask & ~isnan(gt_ori);
			gt_ori(~linemask) = NaN;
			err_i	= ori_error(imgout(linemask),gt_ori(linemask));
			ori_est	= [ori_est; imgout(linemask)];
			ori_gt	= [ori_gt; gt_ori(linemask)];

			% save images to disk
		% 	fname = [outpath,sprintf('%03d_orientation_full.png',i)];
		%		imwrite(uint8(1+255*imgout/pi),cmap,fname);
			fname = [outpath,sprintf('%03d_orientation_masked.png',i)];
				imwrite(uint8(1+255*mod(angle(imgout_masked)*0.5,pi)/pi),ori_cmap,fname);
			fname = [outpath,sprintf('%03d_abs_error.png',i)];
				err_mag = nan(size(linemask)); err_mag(linemask) = abs(err_i);
				imwrite(uint8(1 + 255*err_mag/(pi/2)),mag_cmap,fname);
	end
	
	% update timebar
	timebar(tb,'advance');
end
timebar(tb,'close'); clear tb;

if ~isfield(sampling_args,'prediction_type')
	sampling_args.prediction_type = sampling_args.detection_type;
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
			legend(h,sprintf('%s (%s, %ix%i): mean %0.2f, med %0.2f',...
						strrep(pred_type,'_','\_'),...
						sampling_args.decomp_type,...
						sampling_args.win_size,sampling_args.win_size,...
						es.abs_mean,es.abs_median),...
				'location','southeast');
			title(strrep(pred_job,'_','\_'));
		graph(1); exportfig([outpath,'cumfreq']);

		save([outpath,'errors.mat'],'err','es');
end

filename = [outpath,'error_stats.txt'];
fid = fopen(filename,'w');
	fprintf(fid,'%s\n',evalc('es'));
	fprintf(fid,'%s\n',evalc('time_per_image'));
fclose(fid);

% create a file with the essential properties of the experiment in its name
filename = sprintf('#%s-%s-%ix%i.txt',...
						sampling_args.prediction_type,...
						sampling_args.decomp_type,...
						sampling_args.win_size,sampling_args.win_size);
copyfile([outpath,'error_stats.txt'],[outpath,filename]);

