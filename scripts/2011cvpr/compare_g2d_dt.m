function check_retinogram_errors

clc; clear;

retpath = ...
	[asymmetryroot('shared'),'data\retinograms\DRIVE\test\'];
outpath = ...
	[retpath,'\error_by_property\'];

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

expts = {
	'g2d','analytic';
	'dt','rf_3';
};

for iexpt = 1:length(expts)
	respath = [retpath,'predictions/',...
				expts{iexpt,1},'/',...
				expts{iexpt,2},'/'];

	% get precomputed errors
	load([respath,'errors/orientation_errors.mat']);
	errors_vec = abs(prediction_errs) * 180/pi;
	if strcmp(subset,'centre')
		errors_vec = errors_vec(centres_vec);
	end

	[widths_x(iexpt,:),widths_y(iexpt,:)] = ...
		smoother(widths_vec,errors_vec,100);
	[contrasts_x(iexpt,:),contrasts_y(iexpt,:)] = ...
		smoother(contrasts_vec,errors_vec,100);

	legend_strs{iexpt} = [strrep(expts{iexpt,1},'_','\_'),', ',...
						  strrep(expts{iexpt,2},'_','\_')];
end

graph(1); clf; hold on;
	plot(widths_x',widths_y');
	axis([0,max(widths_vec)*1.05,0,100]);
	set(gca,'box','on');
	xlabel('Line thickness'); ylabel('Abs. angular error (degrees)');
	legend(legend_strs);
	exportfig([outpath,subset,'/thickness_vs_error-summ']);
graph(2); clf; hold on;
	plot(contrasts_x',contrasts_y');
	axis([0,max(contrasts_vec)*1.05,0,100]);
	set(gca,'box','on'); 
	xlabel('Line contrast'); ylabel('Abs. angular error (degrees)');
	legend(legend_strs);
	exportfig([outpath,subset,'/contrast_vs_error-summ']);


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
