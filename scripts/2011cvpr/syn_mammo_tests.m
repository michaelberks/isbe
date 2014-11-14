%**************************************************************************
%***** Script to produce results for CVPR submission Nov 2011 using *******
%****************** the DRIVE retinography dataset ************************
%**************************************************************************
%**************************************************************************
warning('off', 'load_uint8:missing_variables');

clear;

%% 1. Produce response, orientation and scale maps for the analytic methods
retroot = [asymmetryroot('shared'),'data\synthetic_lines\real512\predictions\'];
do_mono =0;
do_g1d = 0;
do_g2d = 1;

if do_mono
    mkdir([retroot 'mono\analytic\orientations']);
    mkdir([retroot 'mono\analytic\responses']);
    mkdir([retroot 'mono\analytic\scales']);
end
if do_g1d
    mkdir([retroot 'g1d\analytic\orientations']);
    mkdir([retroot 'g1d\analytic\responses']);
    mkdir([retroot 'g1d\analytic\scales']);
end
if do_g2d
    mkdir([retroot 'g2d\analytic\orientations']);
    mkdir([retroot 'g2d\analytic\responses']);
    mkdir([retroot 'g2d\analytic\scales']);
end

for ii = 21:100
    %load retinogram and merge RGB channels
    mam = u_load([asymmetryroot('shared'),'data\synthetic_lines\real512\image' zerostr(ii,3) '.mat']);
    
    if do_mono
        %Compute repsonses for mono and save
        [response_map d ori_map scale_map] = monogenic_multiscale(mam, 4, 4, 2, 0.65);
        save_uint8([retroot 'mono\analytic\orientations\' zerostr(ii,3) '_test_ori.mat'], ori_map);
        save_uint8([retroot 'mono\analytic\responses\' zerostr(ii,3) '_test_response.mat'], response_map);
        save_uint8([retroot 'mono\analytic\scales\' zerostr(ii,3) '_test_scale.mat'], scale_map);
    end
    if do_g1d
        %Compute repsonses for mono and save
        [response_map ori_map scale_map] = gaussian_1st_derivative_gradient2(mam, [1 2 4 8]);
        save_uint8([retroot 'g1d\analytic\orientations\' zerostr(ii,3) '_test_ori.mat'], ori_map);
        save_uint8([retroot 'g1d\analytic\responses\' zerostr(ii,3) '_test_response.mat'], response_map);
        save_uint8([retroot 'g1d\analytic\scales\' zerostr(ii,3) '_test_scale.mat'], scale_map);
    end
    if do_g2d
        %Compute repsonses for mono and save
        [response_map ori_map scale_map] = gaussian_clover_line (mam, [1 2 4 8]);
        save_uint8([retroot 'g2d\analytic\orientations\' zerostr(ii,3) '_test_ori.mat'], ori_map);
        save_uint8([retroot 'g2d\analytic\responses\' zerostr(ii,3) '_test_response.mat'], response_map);
        save_uint8([retroot 'g2d\analytic\scales\' zerostr(ii,3) '_test_scale.mat'], scale_map);
    end
    

end

%% 2. Now compute error statistics for all the methods
retroot = [asymmetryroot('shared'),'data\synthetic_lines\real512\'];

%first make sure we've got a complete set of ground truth labels
if 1%~exist([retroot,'orientations\all_gt_orientations.mat'],'file')
	gt_orientations = [];
	line_centres = [];
	pts_per_image = zeros(100,1); 
    for jj = 1:100

		%load line mask and orientation gt
		s = load([retroot 'labels\label' zerostr(jj,3) '.mat']);
		line_mask = s.label == 1;
        centre_mask = s.label_centre & s.label < 2;
        
        line_centres = logical([line_centres; centre_mask(line_mask)]);%#ok

		gt_orientations = [gt_orientations; s.label_orientation(line_mask)]; %#ok
		pts_per_image(jj) = sum(line_mask(:));
    end
    gt_orientations = mod(pi*gt_orientations/180, pi);
	save([retroot,'orientations\all_gt_orientations.mat'], 'gt_orientations', 'pts_per_image', 'line_centres');
end
%%
retroot = [asymmetryroot('shared'),'data\synthetic_lines\real512\'];
load([retroot,'orientations\all_gt_orientations.mat']);

% define jazzy colormaps
ori_cmap = [0 0 0; hsv(255)];
mag_cmap = [0 0 0; jet(255)];        

pred_type{1}  = 'analytic';             pred_decomp{1}  = 'mono';   do_pred(1)  = 0; double_angle(1)  = 0;
pred_type{2}  = 'analytic';             pred_decomp{2}  = 'g1d';    do_pred(2)  = 0; double_angle(2)  = 0;
pred_type{3}  = 'analytic';             pred_decomp{3}  = 'g2d';    do_pred(3)  = 1; double_angle(3)  = 0;
pred_type{4}  = 'rf_3';                 pred_decomp{4}  = 'dt';     do_pred(4)  = 0; double_angle(4)  = 1;
pred_type{5}  = 'rf_3';                 pred_decomp{5}  = 'mono';   do_pred(5)  = 0; double_angle(5)  = 1;
pred_type{6}  = 'rf_3';                 pred_decomp{6}  = 'g1d';    do_pred(6)  = 0; double_angle(6)  = 0;
pred_type{7}  = 'rf_3';                 pred_decomp{7}  = 'g2d';    do_pred(7)  = 0; double_angle(7)  = 1;
pred_type{8}  = 'boosted_regression_3';   pred_decomp{8}  = 'dt';     do_pred(8)  = 0; double_angle(8)  = 0;   
pred_type{9}  = 'boosted_regression_3';   pred_decomp{9}  = 'mono';   do_pred(9)  = 0; double_angle(9)  = 0;
pred_type{10} = 'boosted_regression_3';   pred_decomp{10} = 'g1d';    do_pred(10) = 0; double_angle(10) = 0;
pred_type{11} = 'boosted_regression_3';   pred_decomp{11} = 'g2d';    do_pred(11) = 0; double_angle(11) = 0;
pred_type{12} = 'linear_regression_3';    pred_decomp{12} = 'dt';     do_pred(12) = 0; double_angle(12) = 0;
pred_type{13} = 'linear_regression_3';    pred_decomp{13} = 'mono';   do_pred(13) = 0; double_angle(13) = 0;
pred_type{14} = 'linear_regression_3';    pred_decomp{14} = 'g1d';    do_pred(14) = 0; double_angle(14) = 0;
pred_type{15} = 'linear_regression_3';    pred_decomp{15} = 'g2d';    do_pred(15) = 0; double_angle(15) = 0;

for ii = 1:length(pred_type)
    if ~do_pred(ii); continue; end
    
    err_dir = [retroot 'predictions\' pred_decomp{ii} '\' pred_type{ii} '\errors\'];
    if ~isdir(err_dir); mkdir(err_dir); end
    
    if strcmp(pred_type{ii}, 'analytic')
        ori_dir = [retroot 'predictions\' pred_decomp{ii} '\' pred_type{ii} '\orientations\'];
        res_dir = [retroot 'predictions\' pred_decomp{ii} '\' pred_type{ii} '\responses\'];
        res_list = dir([res_dir '*.mat']);
    else
        ori_dir = [retroot 'predictions\' pred_decomp{ii} '\' pred_type{ii} '\'];
        res_dir = [];
    end
    ori_list = dir([ori_dir '*.mat']);
    
    if length(ori_list) ~= 100
        display([pred_decomp{ii} ' ' pred_type{ii} ': number of ori maps ~= 20; skipping...']);
        continue;
    end
    
    predicted_orientations = [];
    orientation_responses = [];
    for jj = 1:100

        %load line mask and orientation gt
		s = load([retroot 'labels\label' zerostr(jj,3) '.mat']);
		line_mask = s.label == 1;

        %load in orientation map
        ori_map = load_uint8([ori_dir ori_list(jj).name]);
        ori_map(~line_mask) = NaN;

        if isempty(res_dir)
            %compute response from ori_map
            response_map = abs(ori_map);
        else
            %load the response map
            response_map = load_uint8([res_dir res_list(jj).name]);

            %convert the orientation to complex form
            ori_map = exp(complex(0, 2*ori_map));
        end

        if double_angle(ii);
            %double the angle in the complex representation
            ori_map = ori_map.^2;
        end

        %Get the orientation predictions and responses
        predicted_orientations = [predicted_orientations; ori_map(line_mask)]; %#ok
        orientation_responses = [orientation_responses; response_map(line_mask)]; %#ok

    end
    %Now compute the orientation errors for the set
    [prediction_errs all_stats] = ori_error(predicted_orientations, gt_orientations);
    [centre_errs centre_stats] = ori_error(predicted_orientations(line_centres), gt_orientations(line_centres));
    
    %save the predictions and errors
    save([err_dir 'orientation_errors.mat'], 'prediction_errs', 'centre_errs');
    save([err_dir 'orientation_predictions.mat'], 'predicted_orientations');
    save([err_dir 'orientation_responses.mat'], 'orientation_responses');

    %display the results
    display([pred_decomp{ii} ' ' pred_type{ii} ': ' num2str(median(abs(prediction_errs))*180/pi)]);
    
    %produce a txt file with the results
    filename = [err_dir pred_decomp{ii} '_' pred_type{ii} '_error_stats.txt'];
    fid = fopen(filename, 'w');
    fprintf(fid,'%s\n', evalc('all_stats'));
    fprintf(fid,'%s\n', evalc('centre_stats'));
    fclose(fid);
end
%% 3. Draw graphs
retroot = [asymmetryroot('shared'),'data\synthetic_lines\real512\'];
load([retroot,'orientations\all_gt_orientations.mat']);
load([retroot, 'contrasts\centre_gt_contrasts.mat']);

pred_type{1}  = 'analytic';             pred_decomp{1}  = 'mono';   label{1}  = 'Monogenic';
pred_type{2}  = 'analytic';             pred_decomp{2}  = 'g1d';    label{2}  = 'Analytic G''';
pred_type{3}  = 'analytic';             pred_decomp{3}  = 'g2d';    label{3}  = 'Analytic G"';
pred_type{4}  = 'rf_3';                 pred_decomp{4}  = 'dt';     label{4}  = 'Forest DT-CWT';
pred_type{5}  = 'rf_3';                 pred_decomp{5}  = 'mono';   label{5}  = 'Forest Mono.';
pred_type{6}  = 'rf_3';                 pred_decomp{6}  = 'g1d';    label{6}  = 'Forest G''';
pred_type{7}  = 'rf_3';                 pred_decomp{7}  = 'g2d';    label{7}  = 'Forest G"';
pred_type{8}  = 'boosted_regression_3';   pred_decomp{8}  = 'dt';     label{8}  = 'DT-CWT/Boost Reg, 3x3';
pred_type{9}  = 'boosted_regression_3';   pred_decomp{9}  = 'mono';   label{9}  = 'Monogenic/Boost Reg, 3x3';
pred_type{10} = 'boosted_regression_3';   pred_decomp{10} = 'g1d';    label{10} = 'Gauss. 1^{st}/Boost Reg, 3x3';
pred_type{11} = 'boosted_regression_3';   pred_decomp{11} = 'g2d';    label{11} = 'Gauss. 2^{nd}/Roost Reg, 3x3';
pred_type{12} = 'linear_regression_3';    pred_decomp{12} = 'dt';     label{12} = 'LinReg DT-CWT';
pred_type{13} = 'linear_regression_3';    pred_decomp{13} = 'mono';   label{13} = 'Monogenic/Lin Reg, 3x3';
pred_type{14} = 'linear_regression_3';    pred_decomp{14} = 'g1d';    label{14} = 'Gauss. 1^{st}/Lin Reg, 3x3';
pred_type{15} = 'linear_regression_3';    pred_decomp{15} = 'g2d';    label{15} = 'Gauss. 2^{nd}/Lin Reg, 3x3';

% f1 = figure(...
%     'windowstyle', 'normal',...
%     'Units', 'pixels',...
%     'position', [50 50 1000 800],...
%     'PaperPositionMode','auto',...
%     'visible', 'off');
f1 = figure('windowstyle', 'normal','visible', 'off');
graph(f1);
a1 = axes; hold all;
set(a1,'box','on');
ori_leg = cell(0,1);
%title('Synthetic mammogram data'); 
ylabel('Cum. Freq.'); 
xlabel('Orientation error (degrees)');

% title('\fontsize{24} \bf  CDF of orientation estimation errors'); 
% xlabel('\fontsize{20} \bf Orientation error (degrees)');
%ori_leg{1} = '{\fontsize{20} \bf Mean (Median) Angular Errors} \bf \color{white} z_z';
%plot(a1, 0, 0, 'w');

% f2 = figure(...
%     'windowstyle', 'normal',...
%     'Units', 'pixels',...
%     'position', [50 50 1000 800],...
%     'PaperPositionMode','auto',...
%     'visible', 'off');
% a2 = axes; hold all; 
% title({'\fontsize{24} \bf Mean orientation error for the Nth percentile'; '\fontsize{24} \bf of dispersion magnitudes'});
% ylabel('\fontsize{20} \bf Mean orientation error')
% xlabel('\fontsize{20} \bf Percentile of samples sorted by dispersion magnitude');
% pct_leg{1} = '{\fontsize{20} \bf Mean at 50th pcntile} \bf \color{white} z_z';
% plot(a2, 0, 0, 'w');

% f3 = figure(...
%     'windowstyle', 'normal',...
%     'Units', 'pixels',...
%     'position', [50 50 1000 800],...
%     'PaperPositionMode','auto',...
%     'visible', 'off');
% a3 = axes; hold all; 
% title({'\fontsize{24} \bf Median orientation error for the Nth percentile'; '\fontsize{24} \bf of line contrast'});
% ylabel('\fontsize{20} \bf Mean orientation error')
% xlabel('\fontsize{20} \bf Percentile of samples sorted by line contrast');
% con_leg{1} = '{\fontsize{20} \bf Mean at 50th pcntile} \bf \color{white} z_z';
% plot(a3, 0, 0, 'w');

plot_num = 1;
for ii = [3 7 12 4]
    
    err_dir = [retroot 'predictions\' pred_decomp{ii} '\' pred_type{ii} '\errors\'];
    
    %save the predictions and errors
    %load([err_dir 'orientation_predictions.mat'], 'predicted_orientations');
    load([err_dir 'orientation_responses.mat'], 'orientation_responses');
    
    %prediction_errs = abs(angle(predicted_orientations .* conj(gt_orientations))*90/pi);
    
    load([err_dir 'orientation_errors.mat'], 'prediction_errs', 'centre_errs');
    prediction_errs = 180*abs(prediction_errs)/pi;
    centre_errs = 180*abs(centre_errs)/pi;
    
    %get sorted copy of orientation errors
    sorted_ori_errors = sort(prediction_errs);
    
%     %get copy of orientations errors sorted by magnitude of response
%     mag_sorted_errors = sortrows([prediction_errs abs(orientation_responses)], -2);
    
%     %get copy of orientations errors sorted by magnitude of response
%     con_sorted_errors = sortrows([centre_errs line_contrasts_centre], 2);
    
    num_pts = length(prediction_errs);
%     num_pts_c = length(centre_errs);
    ori_cdf = zeros(101,1);
%     mean_pct = zeros(100,1);
%     con_pct = zeros(100,1);
    for jj = 1:100
        x = ceil(jj*num_pts/100);
%         xc = ceil(jj*num_pts_c/100);
        ori_cdf(jj+1) = sorted_ori_errors(x,1);
%         mean_pct(jj) = nanmean(mag_sorted_errors(1:x,1));
%         con_pct(jj) = nanmean(con_sorted_errors(1:xc,1));
    end
    
    plot(a1, ori_cdf, (0:100)/100, 'linewidth', 1);
%     plot(a2, linspace(0,1,100), mean_pct, 'linewidth', 2);
%     plot(a3, linspace(0,1,100), con_pct, 'linewidth', 2);
    
    ori_leg{plot_num} =  label{ii};
    
%     ori_leg{plot_num} = ...
%         ['\fontsize{24} ' label{ii} ': ' sprintf('%3.1f', nanmean(sorted_ori_errors))...
%         ' (' sprintf('%3.1f', naNmedian(sorted_ori_errors)) ') \bf \color{white} z_z']; %#ok
%     pct_leg{plot_num+1} = ...
%         ['\fontsize{20} ' label{ii} ': ' sprintf('%3.1f', mean_pct(50))]; %#ok
%     con_leg{plot_num+1} = ...
%         ['\fontsize{20} ' label{ii} ': ' sprintf('%3.1f', con_pct(50))]; %#ok
    plot_num = plot_num + 1;
end

legend(a1, ori_leg, 'location', 'southeast');
% legend(a2, pct_leg, 'location', 'southeast');
%legend(a3, con_leg, 'location', 'southeast');

set(f1, 'visible', 'on');
% set(f2, 'visible', 'on');
%set(f3, 'visible', 'on');

exportfig([retroot 'figs\cdf_error_mammo.pdf']);
saveas(f1, [retroot 'figs\cdf_error_mammo.fig']);
%print('-dtiff', '-noui', '-painters', f1, '-r300', [retroot 'figs\cdf_error_mammo.pdf']);
% print('-dtiff', '-noui', '-painters', f2, '-r300', [retroot 'figs\response_vs_error_mammo.pdf']);
            
        