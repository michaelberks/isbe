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

% get logical map of pixels with width belwow a threshold
thin_pixels = (widths_vec<=2000);

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
	
for ireg = 1%:length(regressors)
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

		widths_x(ifeat,:) = sort(errors_vec(thin_pixels));
		widths_y(ifeat,:) = (1:sum(thin_pixels))/sum(thin_pixels);
	end

	graph(1); clf; hold on;
		plot(widths_x',widths_y');
		axis([0,max(widths_x(:))*1.05,0,1.02]);
		set(gca,'box','on');
		xlabel('Abs. angular error (degrees)'); ylabel('Cum. Freq.');
		title(strrep(regressors{ireg},'_','\_')); 
		legend(features,'location','southeast');
% 		exportfig([outpath,subset,'/thickness_thresh_2-',(regressors{ireg}]);
% 	graph(2); clf; hold on;
% 		plot(contrasts_x',contrasts_y');
% 		axis([0,max(contrasts_vec)*1.05,0,100]);
% 		set(gca,'box','on'); 
% 		xlabel('Line contrast'); ylabel('Abs. angular error (degrees)');
% 		title(strrep(regressors{ireg},'_','\_')); legend(features);
% 		exportfig([outpath,subset,'/contrast_vs_error-',regressors{ireg}]);
end

