%**************************************************************************
%***** Script to produce results for CVPR submission Nov 2011 using *******
%****************** the DRIVE retinography dataset ************************
%**************************************************************************
%**************************************************************************
warning('off', 'load_uint8:missing_variables');

clear;

%% 1. Produce response, orientation and scale maps for the analytic methods
retroot = [asymmetryroot('shared'),'data\retinograms\DRIVE\test\predictions\orientation\'];
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

for ii = 1:20
    %load retinogram and merge RGB channels
    ret = u_load([asymmetryroot('shared'),'data\retinograms\DRIVE\test\images_extended\' zerostr(ii,2) '_test_ext.mat']);
    %ret = u_load([asymmetryroot('shared'),'data\retinograms\DRIVE\test\images\' zerostr(ii,2) '_test.mat']);
    ret = rgb2gray(ret);
    
    if do_mono
        %Compute repsonses for mono and save
        [response_map d ori_map scale_map] = monogenic_multiscale(ret, 4, 4, 2, 0.65);
        save_uint8([retroot 'mono\analytic\orientations\' zerostr(ii,3) '_test_ori.mat'], ori_map);
        save_uint8([retroot 'mono\analytic\responses\' zerostr(ii,3) '_test_response.mat'], response_map);
        save_uint8([retroot 'mono\analytic\scales\' zerostr(ii,3) '_test_scale.mat'], scale_map);
    end
    if do_g1d
        %Compute repsonses for mono and save
        [response_map ori_map scale_map] = gaussian_1st_derivative_gradient2(ret, [1 2 4 8]);
        save_uint8([retroot 'g1d\analytic\orientations\' zerostr(ii,3) '_test_ori.mat'], ori_map);
        save_uint8([retroot 'g1d\analytic\responses\' zerostr(ii,3) '_test_response.mat'], response_map);
        save_uint8([retroot 'g1d\analytic\scales\' zerostr(ii,3) '_test_scale.mat'], scale_map);
    end
    if do_g2d
        %Compute repsonses for mono and save
        [response_map ori_map scale_map] = gaussian_2nd_derivative_line(ret, [1 2 4 8]);
        save_uint8([retroot 'g2d\analytic\orientations\' zerostr(ii,3) '_test_ori.mat'], ori_map);
        save_uint8([retroot 'g2d\analytic\responses\' zerostr(ii,3) '_test_response.mat'], response_map);
        save_uint8([retroot 'g2d\analytic\scales\' zerostr(ii,3) '_test_scale.mat'], scale_map);
    end
    

end

%% 2. Now compute error statistics for all the methods
retroot = [asymmetryroot('shared'),'data\retinograms\DRIVE\test\'];

%first make sure we've got a complete set of ground truth labels
if ~exist([retroot,'orientations\all_gt_orientations.mat'],'file')
	gt_orientations = [];
	vessel_centres = [];
	pts_per_image = zeros(20,1); 
	for jj = 1:20

		%load foveal mask and vessel mask and indices to vessel centres
		foveal_mask = u_load([retroot 'foveal_masks\' zerostr(jj,2) '_test_f_mask']);
		vessel_mask = u_load([retroot 'vessel_masks\' zerostr(jj,2) '_test_v_mask.mat']);
		vessel_mask(~foveal_mask) = false;

		vessel_centres = logical([vessel_centres; ...
			u_load([retroot 'vessel_masks\centre_idx\' zerostr(jj,2) '_test_centre_idx.mat'])]);%#ok

		%load in the ground truth orientation
		gt_ori = u_load([retroot,'orientations\' zerostr(jj,2) '_ori1.mat']);
		gt_orientations = [gt_orientations; gt_ori(vessel_mask)]; %#ok
		pts_per_image(jj) = sum(vessel_mask(:));
	end
	save([retroot,'orientations\all_gt_orientations.mat'], 'gt_orientations', 'pts_per_image', 'vessel_centres');
end
%%
retroot = [asymmetryroot('shared'),'data\retinograms\DRIVE\test\'];
load([retroot,'orientations\all_gt_orientations.mat']);
pts_im_list = [];
for ii = 1:20
    pts_im_list = [pts_im_list; ii*ones(pts_per_image(ii),1)]; %#ok;
end
pts_im_list(isnan(gt_orientations)) = [];

% define jazzy colormaps
ori_cmap = [0 0 0; hsv(255)];
mag_cmap = [0 0 0; jet(255)];        

pred_type{1}  = 'analytic';             pred_decomp{1}  = 'mono';   do_pred(1)  = 0; double_angle(1)  = 0;
pred_type{2}  = 'analytic';             pred_decomp{2}  = 'g1d';    do_pred(2)  = 0; double_angle(2)  = 0;
pred_type{3}  = 'analytic';             pred_decomp{3}  = 'g2d';    do_pred(3)  = 1; double_angle(3)  = 0;
pred_type{4}  = 'rf_3';                 pred_decomp{4}  = 'dt';     do_pred(4)  = 1; double_angle(4)  = 1;
pred_type{5}  = 'rf_3';                 pred_decomp{5}  = 'mono';   do_pred(5)  = 0; double_angle(5)  = 1;
pred_type{6}  = 'rf_3';                 pred_decomp{6}  = 'g1d';    do_pred(6)  = 0; double_angle(6)  = 0;
pred_type{7}  = 'rf_3';                 pred_decomp{7}  = 'g2d';    do_pred(7)  = 1; double_angle(7)  = 0;
pred_type{8}  = 'boosted_regression_3';   pred_decomp{8}  = 'dt';     do_pred(8)  = 0; double_angle(8)  = 0;   
pred_type{9}  = 'boosted_regression_3';   pred_decomp{9}  = 'mono';   do_pred(9)  = 0; double_angle(9)  = 0;
pred_type{10} = 'boosted_regression_3';   pred_decomp{10} = 'g1d';    do_pred(10) = 0; double_angle(10) = 0;
pred_type{11} = 'boosted_regression_3';   pred_decomp{11} = 'g2d';    do_pred(11) = 0; double_angle(11) = 0;
pred_type{12} = 'linear_regression_3';    pred_decomp{12} = 'dt';     do_pred(12) = 1; double_angle(12) = 0;
pred_type{13} = 'linear_regression_3';    pred_decomp{13} = 'mono';   do_pred(13) = 0; double_angle(13) = 0;
pred_type{14} = 'linear_regression_3';    pred_decomp{14} = 'g1d';    do_pred(14) = 0; double_angle(14) = 0;
pred_type{15} = 'linear_regression_3';    pred_decomp{15} = 'g2d';    do_pred(15) = 1; double_angle(15) = 0;
pred_type{16}  = 'rf_1';                pred_decomp{16}  = 'dt';    do_pred(16) = 0; double_angle(16) = 1;
pred_type{17}  = 'rf_1';                pred_decomp{17}  = 'mono';  do_pred(17) = 0; double_angle(17) = 0;
pred_type{18}  = 'rf_1';                pred_decomp{18}  = 'g1d';   do_pred(18) = 0; double_angle(18) = 0;
pred_type{19}  = 'rf_1';                pred_decomp{19}  = 'g2d';   do_pred(19) = 0; double_angle(19) = 0;
pred_type{end+1}  = 'boosted_regression_1';   pred_decomp{end+1}  = 'dt';     do_pred(end+1)  = 0; double_angle(end+1)  = 0;   
pred_type{end+1}  = 'boosted_regression_1';   pred_decomp{end+1}  = 'mono';   do_pred(end+1)  = 0; double_angle(end+1)  = 0;
pred_type{end+1} = 'boosted_regression_1';   pred_decomp{end+1} = 'g1d';    do_pred(end+1) = 0; double_angle(end+1) = 0;
pred_type{end+1} = 'boosted_regression_1';   pred_decomp{end+1} = 'g2d';    do_pred(end+1) = 0; double_angle(end+1) = 0;
pred_type{end+1} = 'linear_regression_1';    pred_decomp{end+1} = 'dt';     do_pred(end+1) = 0; double_angle(end+1) = 0;
pred_type{end+1} = 'linear_regression_1';    pred_decomp{end+1} = 'mono';   do_pred(end+1) = 0; double_angle(end+1) = 0;
pred_type{end+1} = 'linear_regression_1';    pred_decomp{end+1} = 'g1d';    do_pred(end+1) = 0; double_angle(end+1) = 0;
pred_type{end+1} = 'linear_regression_1';    pred_decomp{end+1} = 'g2d';    do_pred(end+1) = 0; double_angle(end+1) = 0;
pred_type{end+1} = 'rf_3';    pred_decomp{end+1} = 'dtg2';    do_pred(end+1) = 1; double_angle(end+1) = 0;
pred_type{end+1} = 'rf_3';    pred_decomp{end+1} = 'g12d';    do_pred(end+1) = 1; double_angle(end+1) = 0;
pred_type{end+1} = 'rf_3c';    pred_decomp{end+1} = 'dt';    do_pred(end+1) = 1; double_angle(end+1) = 0;
pred_type{end+1} = 'analytic';    pred_decomp{end+1} = 'post_classc_g2d';    do_pred(end+1) = 1; double_angle(end+1) = 0;
pred_type{end+1} = 'rf_3';    pred_decomp{end+1} = 'jim';    do_pred(end+1) = 1; double_angle(end+1) = 0;

for ii = 3%length(pred_type)%[3 4 7 12 15]%%1:
    if ~do_pred(ii); continue; end
    
    err_dir = [retroot 'predictions\orientation\' pred_decomp{ii} '\' pred_type{ii} '\errors\'];
    if ~isdir(err_dir); mkdir(err_dir); end
    
    if strcmp(pred_type{ii}, 'analytic')
        ori_dir = [retroot 'predictions\orientation\' pred_decomp{ii} '\' pred_type{ii} '\orientations\'];
        res_dir = [retroot 'predictions\orientation\' pred_decomp{ii} '\' pred_type{ii} '\responses\'];
        res_list = dir([res_dir '*.mat']);
    else
        ori_dir = [retroot 'predictions\orientation\' pred_decomp{ii} '\' pred_type{ii} '\'];
        res_dir = [];
    end
    ori_list = dir([ori_dir '*.mat']);
    
    if length(ori_list) ~= 20
        display([pred_decomp{ii} ' ' pred_type{ii} ': number of ori maps ~= 20; skipping...']);
        continue;
    end
    
    predicted_orientations = [];
    orientation_responses = [];
    for jj = 1:20

        %load foveal mask and vessel mask and indices to vessel centres
        foveal_mask = u_load([retroot 'foveal_masks\' zerostr(jj,2) '_test_f_mask']);
        vessel_mask = u_load([retroot 'vessel_masks\' zerostr(jj,2) '_test_v_mask.mat']);
        vessel_mask(~foveal_mask) = false;

        %load in orientation map
        ori_map = load_uint8([ori_dir ori_list(jj).name]);
        ori_map(~vessel_mask) = NaN;

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
        predicted_orientations = [predicted_orientations; ori_map(vessel_mask)]; %#ok
        orientation_responses = [orientation_responses; response_map(vessel_mask)]; %#ok
        
        %Save images of the predictions and the error magnitudes
        error_map = nan(size(vessel_mask));
        error_map(vessel_mask) = angle(ori_map(vessel_mask) .* ...
            conj(gt_orientations(sum(pts_per_image(1:jj-1))+1:sum(pts_per_image(1:jj)))))/2;
        
        save([err_dir zerostr(jj,2) '_error_map.mat'], 'error_map');
		imwrite(uint8(1+255*mod(angle(ori_map)/2,pi)/pi), ori_cmap,...
            [err_dir zerostr(jj,2) '_orientation_masked.png']);
        imwrite(uint8(1 + 255*abs(error_map)/(pi/2)), mag_cmap,...
            [err_dir zerostr(jj,2) '_abs_error.png']);

    end
    %Now compute the orientation errors for the set
    [prediction_errs all_stats] = ori_error(predicted_orientations, gt_orientations);
    [centre_errs centre_stats] = ori_error(predicted_orientations(vessel_centres), gt_orientations(vessel_centres));
    
    %Compute bootstrap errors
    median_error_bt = zeros(2000,1);
    for jj = 1:2000
        boot_samp = ceil(20*rand(20,1));
        pts_idx = ismember(pts_im_list, boot_samp);
        median_error_bt(jj) = median(abs(prediction_errs(pts_idx)));
    end
    median_ci = 180*grpstats(median_error_bt, ones(2000,1), 'meanci')/pi;

    %save the predictions and errors
    save([err_dir 'orientation_errors.mat'], 'prediction_errs', 'centre_errs');
    save([err_dir 'orientation_predictions.mat'], 'predicted_orientations');
    save([err_dir 'orientation_responses.mat'], 'orientation_responses');

    %display the results
    display([pred_decomp{ii} ' ' pred_type{ii} ': ' num2str(median(abs(prediction_errs))*180/pi) ...
        ' (' num2str(median_ci(1)) ', ' num2str(median_ci(2)) ')' ]);
    
    %produce a txt file with the results
    filename = [err_dir pred_decomp{ii} '_' pred_type{ii} '_error_stats.txt'];
    fid = fopen(filename, 'w');
    fprintf(fid,'%s\n', evalc('all_stats'));
    fprintf(fid,'%s\n', evalc('centre_stats'));
    fclose(fid);
end
%% 3. Draw graphs
retroot = [asymmetryroot('shared'),'data\retinograms\DRIVE\test\'];
load([retroot,'orientations\all_gt_orientations.mat']);

pred_type{1}  = 'analytic';             pred_decomp{1}  = 'mono';   label{1}  = 'Monogenic';
pred_type{2}  = 'analytic';             pred_decomp{2}  = 'g1d';    label{2}  = 'Gauss. 1^{st} deriv';
pred_type{3}  = 'analytic';             pred_decomp{3}  = 'g2d';    label{3}  = 'Gauss. 2^{nd} deriv';
pred_type{4}  = 'rf_3';                 pred_decomp{4}  = 'dt';     label{4}  = 'DT-CWT/RF, 3x3';
pred_type{5}  = 'rf_3';                 pred_decomp{5}  = 'mono';   label{5}  = 'Monogenic/RF, 3x3';
pred_type{6}  = 'rf_3';                 pred_decomp{6}  = 'g1d';    label{6}  = 'Gauss. 1^{st}/RF, 3x3';
pred_type{7}  = 'rf_3';                 pred_decomp{7}  = 'g2d';    label{7}  = 'Gauss. 2^{nd}/RF, 3x3';
pred_type{8}  = 'boosted_regression';   pred_decomp{8}  = 'dt';     label{8}  = 'DT-CWT/Boost Reg, 3x3';
pred_type{9}  = 'boosted_regression';   pred_decomp{9}  = 'mono';   label{9}  = 'Monogenic/Boost Reg, 3x3';
pred_type{10} = 'boosted_regression';   pred_decomp{10} = 'g1d';    label{10} = 'Gauss. 1^{st}/Boost Reg, 3x3';
pred_type{11} = 'boosted_regression';   pred_decomp{11} = 'g2d';    label{11} = 'Gauss. 2^{nd}/Roost Reg, 3x3';
pred_type{12} = 'linear_regression';    pred_decomp{12} = 'dt';     label{12} = 'DT-CWT/Lin Reg, 3x3';
pred_type{13} = 'linear_regression';    pred_decomp{13} = 'mono';   label{13} = 'Monogenic/Lin Reg, 3x3';
pred_type{14} = 'linear_regression';    pred_decomp{14} = 'g1d';    label{14} = 'Gauss. 1^{st}/Lin Reg, 3x3';
pred_type{15} = 'linear_regression';    pred_decomp{15} = 'g2d';    label{15} = 'Gauss. 2^{nd}/Lin Reg, 3x3';
pred_type{16}  = 'rf_1';                pred_decomp{16}  = 'dt';    label{16} = 'DT-CWT/RF, 1x1';
pred_type{17}  = 'rf_1';                pred_decomp{17}  = 'mono';  label{17} = 'Monogenic/RF, 1x1';
pred_type{18}  = 'rf_1';                pred_decomp{18}  = 'g1d';   label{18} = 'Gauss. 1^{st}/RF, 1x1';
pred_type{19}  = 'rf_1';                pred_decomp{19}  = 'g2d';   label{19} = 'Gauss. 2^{nd}/RF, 1x1';
pred_type{20}  = 'rf_3';                pred_decomp{20}  = 'dtg2';  label{20} = 'DT-CWT+Gauss. 2^{nd}/RF, 3x3';
pred_type{21} =  'rf_3';                pred_decomp{21} =  'g12d';  label{21} = 'Gauss. 1^{st} + 2^{nd}/RF, 3x3';

f1 = figure(...
    'windowstyle', 'normal',...
    'Units', 'pixels',...
    'position', [50 50 1000 800],...
    'PaperPositionMode','auto',...
    'visible', 'off');
a1 = axes; hold all; 
title('\fontsize{24} \bf  CDF of orientation estimation errors'); 
xlabel('\fontsize{20} \bf Orientation error (degrees)');
%ori_leg{1} = '{\fontsize{20} \bf Mean (Median) Angular Errors} \bf \color{white} z_z';
%plot(a1, 0, 0, 'w');

f2 = figure(...
    'windowstyle', 'normal',...
    'Units', 'pixels',...
    'position', [50 50 1000 800],...
    'PaperPositionMode','auto',...
    'visible', 'off');
a2 = axes; hold all; 
title({'\fontsize{24} \bf Mean orientation error for the Nth percentile'; '\fontsize{24} \bf of dispersion magnitudes'});
ylabel('\fontsize{20} \bf Mean orientation error')
xlabel('\fontsize{20} \bf Percentile of samples sorted by dispersion magnitude');
%pct_leg{1} = '{\fontsize{20} \bf Mean at 50th pcntile} \bf \color{white} z_z';
%plot(a2, 0, 0, 'w');

plot_num = 1;
for ii = [4 7 5 6 3 21]
    
    err_dir = [retroot 'predictions\orientation\' pred_decomp{ii} '\' pred_type{ii} '\errors\'];
    
    %save the predictions and errors
    load([err_dir 'orientation_predictions.mat'], 'predicted_orientations');
    load([err_dir 'orientation_responses.mat'], 'orientation_responses');
    
    prediction_errs = abs(angle(predicted_orientations .* conj(gt_orientations))*90/pi);
    
    %get sorted copy of orientation errors
    sorted_ori_errors = sort(prediction_errs);
    
    %get copy of orientations errors sorted by magnitude of response
    mag_sorted_errors = sortrows(abs([prediction_errs orientation_responses]), -2);
    
    num_pts = length(prediction_errs);
    ori_cdf = zeros(101,1);
    mean_pct = zeros(100,1);
    for jj = 1:100
        x = ceil(jj*num_pts/100);
        ori_cdf(jj+1) = sorted_ori_errors(x,1);
        mean_pct(jj) = naNmedian(mag_sorted_errors(1:x,1));
    end
    
    plot(a1, ori_cdf, (0:100)/100, 'linewidth', 2);
    plot(a2, linspace(0,1,100), mean_pct, 'linewidth', 2);
    
    ori_leg{plot_num} = ['\fontsize{24} ' label{ii}]; %#ok
    pct_leg{plot_num} = ['\fontsize{20} ' label{ii}]; %#ok
    plot_num = plot_num + 1;
end

legend(a1, ori_leg, 'location', 'southeast');
legend(a2, pct_leg, 'location', 'southeast');

set(f1, 'visible', 'on');
set(f2, 'visible', 'on');
%%
retroot = [asymmetryroot('shared'),'data\retinograms\DRIVE\test\'];
load([retroot,'retinogram_properties.mat']);

pred_type{1}  = 'analytic';             pred_decomp{1}  = 'mono';   label{1}  = 'Monogenic';
pred_type{2}  = 'analytic';             pred_decomp{2}  = 'g1d';    label{2}  = 'Gauss. 1^{st} deriv';
pred_type{3}  = 'analytic';             pred_decomp{3}  = 'g2d';    label{3}  = 'Gauss. 2^{nd} deriv';
pred_type{4}  = 'rf_3';                 pred_decomp{4}  = 'dt';     label{4}  = 'DT-CWT/RF, 3x3';
pred_type{5}  = 'rf_3';                 pred_decomp{5}  = 'mono';   label{5}  = 'Monogenic/RF, 3x3';
pred_type{6}  = 'rf_3';                 pred_decomp{6}  = 'g1d';    label{6}  = 'Gauss. 1^{st}/RF, 3x3';
pred_type{7}  = 'rf_3';                 pred_decomp{7}  = 'g2d';    label{7}  = 'Gauss. 2^{nd}/RF, 3x3';
pred_type{8}  = 'boosted_regression';   pred_decomp{8}  = 'dt';     label{8}  = 'DT-CWT/Boost Reg, 3x3';
pred_type{9}  = 'boosted_regression';   pred_decomp{9}  = 'mono';   label{9}  = 'Monogenic/Boost Reg, 3x3';
pred_type{10} = 'boosted_regression';   pred_decomp{10} = 'g1d';    label{10} = 'Gauss. 1^{st}/Boost Reg, 3x3';
pred_type{11} = 'boosted_regression';   pred_decomp{11} = 'g2d';    label{11} = 'Gauss. 2^{nd}/Roost Reg, 3x3';
pred_type{12} = 'linear_regression';    pred_decomp{12} = 'dt';     label{12} = 'DT-CWT/Lin Reg, 3x3';
pred_type{13} = 'linear_regression';    pred_decomp{13} = 'mono';   label{13} = 'Monogenic/Lin Reg, 3x3';
pred_type{14} = 'linear_regression';    pred_decomp{14} = 'g1d';    label{14} = 'Gauss. 1^{st}/Lin Reg, 3x3';
pred_type{15} = 'linear_regression';    pred_decomp{15} = 'g2d';    label{15} = 'Gauss. 2^{nd}/Lin Reg, 3x3';
pred_type{16}  = 'rf_1';                pred_decomp{16}  = 'dt';    label{16} = 'DT-CWT/RF, 1x1';
pred_type{17}  = 'rf_1';                pred_decomp{17}  = 'mono';  label{17} = 'Monogenic/RF, 1x1';
pred_type{18}  = 'rf_1';                pred_decomp{18}  = 'g1d';   label{18} = 'Gauss. 1^{st}/RF, 1x1';
pred_type{19}  = 'rf_1';                pred_decomp{19}  = 'g2d';   label{19} = 'Gauss. 2^{nd}/RF, 1x1';
pred_type{20}  = 'rf_3';                pred_decomp{20}  = 'dtg2';  label{20} = 'DT-CWT+Gauss. 2^{nd}/RF, 3x3';
pred_type{21} =  'rf_3';                pred_decomp{21} =  'g12d';  label{21} = 'Gauss. 1^{st} + 2^{nd}/RF, 3x3';

vessel_widths = [];
vessel_centres = [];
for ii = 1:20
    vessel_widths = [vessel_widths; line_widths{ii}];
    vessel_centres = logical([vessel_centres; centre_inds{ii}]);
end

f1 = figure(...
    'windowstyle', 'normal',...
    'Units', 'pixels',...
    'position', [50 50 1000 800],...
    'PaperPositionMode','auto',...
    'visible', 'off');
a1 = axes; hold all; 
title('\fontsize{24} \bf  CDF of orientation estimation errors at vessel centres'); 
xlabel('\fontsize{20} \bf Orientation error (degrees)');
ori_leg_c = cell(0,1);
%ori_leg{1} = '{\fontsize{20} \bf Mean (Median) Angular Errors} \bf \color{white} z_z';
%plot(a1, 0, 0, 'w');

f2 = figure(...
    'windowstyle', 'normal',...
    'Units', 'pixels',...
    'position', [50 50 1000 800],...
    'PaperPositionMode','auto',...
    'visible', 'off');
a2 = axes; hold all; 
title({'\fontsize{24} \bf CDF of orientation estimation errors at vessel centres'; '\fontsize{24} \bf (thin vessels)'});
xlabel('\fontsize{20} \bf Orientation error (degrees)');
ori_leg_t = cell(0,1);
%pct_leg{1} = '{\fontsize{20} \bf Mean at 50th pcntile} \bf \color{white} z_z';
%plot(a2, 0, 0, 'w');

plot_num = 1;
for ii = [4 7 5 6 3 21]
    
    err_dir = [retroot 'predictions\orientation\' pred_decomp{ii} '\' pred_type{ii} '\errors\'];
    
    %load the predictions and errors
    load([err_dir 'orientation_errors.mat']);
    
    centre_errs = prediction_errs(vessel_centres);
    centre_widths = vessel_widths(vessel_centres);
    centre_errs = abs(centre_errs)*180/pi;
    thin_vessels = centre_widths <= 2;
    thin_errs = centre_errs(thin_vessels);
    
    %get sorted copy of orientation errors
    sorted_errors_c = sort(centre_errs);
    sorted_errors_t = sort(thin_errs);
    
    num_pts_c = length(centre_errs);
    num_pts_t = length(thin_errs);
    ori_cdf_c = zeros(101,1);
    ori_cdf_t = zeros(101,1);
    for jj = 1:100
        x_c = ceil(jj*num_pts_c/100);
        x_t = ceil(jj*num_pts_t/100);
        ori_cdf_c(jj+1) = sorted_errors_c(x_c,1);
        ori_cdf_t(jj+1) = sorted_errors_t(x_t,1);
    end
    
    plot(a1, ori_cdf_c, (0:100)/100, 'linewidth', 2);
    plot(a2, ori_cdf_t, (0:100)/100, 'linewidth', 2);
    
    ori_leg_c{plot_num} = ['\fontsize{24} ' label{ii}]; %#ok
    ori_leg_t{plot_num} = ['\fontsize{24} ' label{ii}]; %#ok
    plot_num = plot_num + 1;
end

legend(a1, ori_leg_c, 'location', 'southeast');
legend(a2, ori_leg_t, 'location', 'southeast');

set(f1, 'visible', 'on');
set(f2, 'visible', 'on');
            
        