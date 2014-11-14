% Script for predicting orientation in synthetic line images and producing
% stuff for the BMVC paper.
%
% Note that for these experiments I'm using regressors trained with 20k
% data points, usually on a 1x1 grid though some will use 3x3

close all hidden;
clear variables;

%% define predictor type and load
results_ind = '';
% pred_type	= 'linear_regression'; pred_job	= '313751/'; % 200k dt, generated
pred_type	= 'orientation'; pred_job = '286712/'; % 200k dt, generated


% default sampling arguments
sampling_args = struct(...
	'decomp_type','',...
	'rgb_channel','',...
	'win_size',-1,...
	'prediction_type',pred_type);

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


%% process test images

% define and create output folder
retpath = 'mammograms/2004_screening_processed/mass_roi';
retroot = [asymmetryroot,'data/',retpath,'/'];
imgpath = [retpath];
imgroot = [asymmetryroot,'data/',imgpath,'/'];
outpath = strrep([imgroot,'results/',pred_job,'/',results_ind],'\','/');
outpath = strrep(outpath,'//','/');
if ~exist(outpath,'dir'), mkdir(outpath); end

if exist(pred_path,'dir')
	% copy parameter list to results dir
	if exist([pred_path,'args.txt'],'file')
		copyfile([pred_path,'args.txt'],outpath);
	end
end

% imgdir = dir([imgroot,'*.mat']);
% for i = 1:length(imgdir)
% 	img = load_uint8([imgroot,imgdir(i).name]);
% 	save([imgroot,'double/',imgdir(i).name],'img');
% end
imgpath = [imgpath,'/double'];
imgroot = [asymmetryroot,'data/',imgpath,'/'];
outpath = strrep([imgroot,'results/',pred_job,'/',results_ind],'\','/');
outpath = strrep(outpath,'//','/');
if ~exist(outpath,'dir'), mkdir(outpath); end

view = '024RCC';
% view = [];

% resdir = dir([outpath,'*_class.mat']);
% if length(resdir)<2
% 	% classify set of test images
% 	tic;
% 	classify_image_set(pred_job,imgpath,...
% 		'view',view,...
% 		'forest_dir', sprintf('models/mammograms/%s',pred_type), ...
% 		'use_nag', ~ispc);
% 	T = toc;
% end

%% generate PNG images from saved predictions

% define jazzy colormaps
ori_cmap = [0 0 0; hsv(255)];
mag_cmap = [0 0 0; jet(255)];

resdir = dir([outpath,'*_class.mat']);
tb = timebar('limit',length(resdir),'title','Converting images');
for i = 1:length(resdir)
	% load predicted output back in
	imgout = load_uint8([outpath,resdir(i).name]);

	[p,outbase,e] = fileparts(resdir(i).name);
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

		case {'orientation','linear_regression','logistic_regression','boosted_regression'},
			% angle needs doubling for 'orientation'
			if strcmp(pred_type,'orientation')
% 				imgout = abs(imgout).^2 .* exp(complex(0,angle(imgout)*2));
% 				imgout = complex2rgb(imgout.^2);
				imgout2(:,:,1) = 
			end				
			
			% save images to disk
		% 	fname = [outpath,sprintf('%03d_orientation_full.png',i)];
		%		imwrite(uint8(1+255*imgout/pi),cmap,fname);
			fname = [outpath,sprintf('%s_ori_masked.png',outbase)];
			if size(imgout,3)==1
				imwrite(uint8(1+255*mod(angle(imgout)/2,pi)/pi),ori_cmap,fname);
			else
				imwrite(uint8(255*imgout),fname);
			end
	end
	
	% update timebar
	timebar(tb,'advance');
end
timebar(tb,'close'); clear tb;
