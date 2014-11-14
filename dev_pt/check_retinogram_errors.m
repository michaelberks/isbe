function check_retinogram_errors

clc; clear;

retpath = ...
	[asymmetryroot('shared'),'data\retinograms\DRIVE\test\'];
outpath = ...
	'S:\projects\mammography\matlab\papers\2011cvpr\figs\retina\error_by_property\';

% load line properties and pool measurements
load([retpath,'retinogram_properties.mat']);
widths_vec = cat(1,line_widths{:});
contrasts_vec = cat(1,line_contrasts{:});
centres_vec = cat(1,centre_inds{:});

% consider only centrelines
subset = 'vessel';
if 1
	widths_vec = widths_vec(centres_vec);
	contrasts_vec = contrasts_vec(centres_vec);
	subset = 'centre';
end

features = {
	'g1d';
	'g2d';
	'mono';
	'dt';
};

regressors = {
	'linear_regression_1';
	'boosted_regression_1';
	'rf_1';
};
	
for ireg = 1:length(regressors)
	for ifeat = 1:length(features)
		respath = [retpath,'predictions/',...
					features{ifeat},'/',...
					regressors{ireg},'/'];

		% get precomputed errors
		load([respath,'errors/orientation_errors.mat']);
		errors_vec = abs(prediction_errs) * 180/pi;
		if strcmp(subset,'centre')
			errors_vec = errors_vec(centres_vec);
		end

		[widths_x(ifeat,:),widths_y(ifeat,:)] = ...
			smoother(widths_vec,errors_vec,100);
		[contrasts_x(ifeat,:),contrasts_y(ifeat,:)] = ...
			smoother(contrasts_vec,errors_vec,100);
	end

	graph(1); clf; hold on;
		plot(widths_x',widths_y');
		axis([0,max(widths_vec)*1.05,0,100]);
		set(gca,'box','on');
		xlabel('Line thickness'); ylabel('Abs. angular error (degrees)');
		title(strrep(regressors{ireg},'_','\_')); legend(features);
		exportfig([outpath,subset,'/thickness_vs_error-',regressors{ireg}]);
	graph(2); clf; hold on;
		plot(contrasts_x',contrasts_y');
		axis([0,max(contrasts_vec)*1.05,0,100]);
		set(gca,'box','on'); 
		xlabel('Line contrast'); ylabel('Abs. angular error (degrees)');
		title(strrep(regressors{ireg},'_','\_')); legend(features);
		exportfig([outpath,subset,'/contrast_vs_error-',regressors{ireg}]);
end


function [centres,smoothed_y] = smoother(x,y,N,kernel,varargin)

if nargin<3, N = 100; end
if nargin<4, kernel = 'gaussian'; end

[sorted_x,inds] = sort(x);
sorted_y = y(inds);

centres = linspace(sorted_x(1),sorted_x(end),N);
smoothed_y = zeros(size(centres));

switch kernel
	case 'gaussian',
		% gaussian kernel density estimate
		
		if length(varargin)<1, w_sigma_scale = 3; end
		
		w_sigma = (centres(2)-centres(1)) * w_sigma_scale;
		for i = 1:N
			normed = (sorted_x - centres(i)) / w_sigma;
			weights = exp(-0.5*normed.*normed);
			smoothed_y(i) = sum(weights.*sorted_y) / sum(weights);
		end

	case 'cummean',
		% cumulative mean
		for i = 1:N
			weights = (sorted_x <= centres(i));
			smoothed_y(i) = sum(weights.*sorted_y) / sum(weights);
		end
		
	otherwise,
		error(['Unknown kernel: ',kernel]);
end
