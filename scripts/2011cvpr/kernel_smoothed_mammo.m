function kernel_smoothed_mammo

%clc; clear;

retroot = [asymmetryroot('shared'),'data\retinograms\DRIVE\test\'];
load([retroot,'orientations\all_gt_orientations.mat']);
%load([retroot, 'contrasts\centre_gt_contrasts.mat']);

pred_type{1}  = 'analytic';             pred_decomp{1}  = 'mono';   label{1}  = 'Monogenic';
pred_type{2}  = 'analytic';             pred_decomp{2}  = 'g1d';    label{2}  = 'Gauss. 1^{st} deriv';
pred_type{3}  = 'analytic';             pred_decomp{3}  = 'g2d';    label{3}  = 'Gauss. 2^{nd} deriv';
pred_type{4}  = 'rf_3';                 pred_decomp{4}  = 'dt';     label{4}  = 'Forest DT-CWT';
pred_type{5}  = 'rf_3';                 pred_decomp{5}  = 'mono';   label{5}  = 'Forest Mono.';
pred_type{6}  = 'rf_3';                 pred_decomp{6}  = 'g1d';    label{6}  = 'Forest G''';
pred_type{7}  = 'rf_3';                 pred_decomp{7}  = 'g2d';    label{7}  = 'Forest G"';
pred_type{8}  = 'boosted_regression_3';   pred_decomp{8}  = 'dt';     label{8}  = 'DT-CWT/Boost Reg, 3x3';
pred_type{9}  = 'boosted_regression_3';   pred_decomp{9}  = 'mono';   label{9}  = 'Monogenic/Boost Reg, 3x3';
pred_type{10} = 'boosted_regression_3';   pred_decomp{10} = 'g1d';    label{10} = 'Gauss. 1^{st}/Boost Reg, 3x3';
pred_type{11} = 'boosted_regression_3';   pred_decomp{11} = 'g2d';    label{11} = 'Gauss. 2^{nd}/Roost Reg, 3x3';
pred_type{12} = 'linear_regression_3';    pred_decomp{12} = 'dt';     label{12} = 'DT-CWT/Lin Reg, 3x3';
pred_type{13} = 'linear_regression_3';    pred_decomp{13} = 'mono';   label{13} = 'Monogenic/Lin Reg, 3x3';
pred_type{14} = 'linear_regression_3';    pred_decomp{14} = 'g1d';    label{14} = 'Gauss. 1^{st}/Lin Reg, 3x3';
pred_type{15} = 'linear_regression_3';    pred_decomp{15} = 'g2d';    label{15} = 'Gauss. 2^{nd}/Lin Reg, 3x3';

% plot_num = 1;
% contrasts_x = [];
% contrasts_y = [];
% for ii = [4 7 12 3]
%     
%     err_dir = [retroot 'predictions\' pred_decomp{ii} '\' pred_type{ii} '\errors\'];
%     
%     %save the predictions and errors
%     load([err_dir 'orientation_errors.mat'], 'centre_errs');
%     centre_errs = 180*abs(centre_errs)/pi;
%     
%     
%     [contrasts_x(plot_num,:), contrasts_y(plot_num,:)] = ...
%         smoother(line_contrasts_centre, centre_errs, 100);
%         
%     plot_num = plot_num + 1;
% 
%     
%     
% end
% 
% graph(1); clf; hold on;
% plot(contrasts_x',contrasts_y');
% axis([0,max(line_contrasts_centre)*1.05,0,100]);
% set(gca,'box','on'); 
% xlabel('Line contrast'); ylabel('Abs. angular error (degrees)');
% legend(label([4 7 12 3]));
% exportfig([retroot 'figs\contrast_vs_error_mammo.pdf']);

plot_num = 1;
responses_x = [];
responses_y = [];
for ii = [3 6 7 5 4]
    
    err_dir = [retroot 'predictions\orientation\' pred_decomp{ii} '\' pred_type{ii} '\errors\'];
    
    %save the predictions and errors
    %load([err_dir 'orientation_predictions.mat'], 'predicted_orientations');
    load([err_dir 'orientation_responses.mat'], 'orientation_responses');
    load([err_dir 'orientation_errors.mat'], 'prediction_errs');
    
    orientation_responses(isnan(gt_orientations)) = [];
    prediction_errs = 180*abs(prediction_errs)/pi;

    if strcmpi(pred_type{ii}, 'analytic');
        orientation_responses = abs(orientation_responses);
        orientation_responses = orientation_responses / max(orientation_responses);
    end
    
    [responses_x(plot_num,:), responses_y(plot_num,:)] = ...
        smoother(orientation_responses, prediction_errs, 100);
        
    plot_num = plot_num + 1;
  
end

f1 = figure('windowstyle', 'normal');
graph(f1); hold on;
plot(responses_x',responses_y');
axis([0,1,0,45]);
set(gca,'box','on'); 
xlabel('Magnitude of predicted orientation vector'); ylabel('Abs. angular error (degrees)');
legend(label([6 7 5 4]), 'location', 'southwest');
exportfig([retroot 'figs\response_vs_error_mammo.pdf']);

%save([retroot 'figs\response_vs_error_mammo.mat'], 'responses_x', 'responses_y');


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
