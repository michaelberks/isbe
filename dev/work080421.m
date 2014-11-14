%work 21/04/2008
%%
% Experiments looking at the difference in pyramid coefficients at CLS
% pixels depending on the orientation of the CLS and sub-band
%
% For any sub-band we can measure the average coefficient value for:
% 1) all pixels
% 2) all CLS pixels
% 3) all non-CLS pixels
% 4) all aligned CLS pixels
% 5) all non-aligned CLS pixels
% 6) all non-aligned or non-CLS pixels

% pyr_all = zeros(53,5,3);
% pyr_all_cls = zeros(53,5,3);
% pyr_non_cls = zeros(53,5,3);
% pyr_aligned_cls = zeros(53,5,3);
% pyr_nonaligned_cls = zeros(53,5,3);
% pyr_all_nonaligned_cls = zeros(53,5,3);
% 
% dt_all = zeros(53,6,3);
% dt_all_cls = zeros(53,6,3);
% dt_non_cls = zeros(53,6,3);
% dt_aligned_cls = zeros(53,6,3);
% dt_nonaligned_cls = zeros(53,6,3);
% dt_all_nonaligned_cls = zeros(53,6,3);

% dt_abs_all = zeros(53,6,3);
% dt_abs_all_cls = zeros(53,6,3);
% dt_abs_non_cls = zeros(53,6,3);
% dt_abs_aligned_cls = zeros(53,6,3);
% dt_abs_nonaligned_cls = zeros(53,6,3);
% dt_abs_all_nonaligned_cls = zeros(53,6,3);

pyr_ttest_ab = zeros(53,5,3);
pyr_ttest_bc = zeros(53,5,3);
pyr_ttest_ac = zeros(53,5,3);
%
dt_ttest_real_ab = zeros(53,5,3);
dt_ttest_real_bc = zeros(53,5,3);
dt_ttest_real_ac = zeros(53,5,3);
dt_ttest_imag_ab = zeros(53,5,3);
dt_ttest_imag_bc = zeros(53,5,3);
dt_ttest_imag_ac = zeros(53,5,3);

orientations = [Inf -Inf; 0 pi; pi 2*pi; 2*pi 3*pi; -3*pi -2*pi; -2*pi -pi; -pi 0] / 6;

for ii = 1:53 %53 = number of masses
    %Load CLS result for each mass
    cls = u_load(['C:\isbe\dev\background\cls\mass_2\mass_cls_', zerostr(ii,3), '.mat']);
    
    %Generate oriented CLS map
    [cls_map] =  mb_cls_map_orientations(cls);
    
    %Load pyramid
    pyr = u_load(['C:\isbe\dev\background\pyramid5\mass_2\mass', zerostr(ii,3) ,'_pyramid.mat']);
    
    %Compute values for each pixel typ in each sub-band
    for lev = 1:3
        for ori = 1:5           
            
            temp = (pyr{lev+1,ori});
            %Compute means for each type
%             pyr_all(ii, ori, lev) = mean(temp(:));
%             pyr_all_cls(ii, ori, lev) = mean(temp(cls_map{lev}>0));
%             pyr_non_cls(ii, ori, lev) = mean(temp(~cls_map{lev}));
%             pyr_aligned_cls(ii, ori, lev) = mean(temp(cls_map{lev} == ori));
%             pyr_nonaligned_cls(ii, ori, lev) = mean(temp(cls_map{lev} &  cls_map{lev} ~= ori));
%             pyr_all_nonaligned_cls(ii, ori, lev) = mean(temp(cls_map{lev} ~= ori));
            
            [h p] = ttest2(temp(~cls_map{lev}), temp(cls_map{lev} &  cls_map{lev} ~= ori));
            pyr_ttest_ab(ii, ori, lev) = p;
            [h p] = ttest2(temp(cls_map{lev} == ori), temp(cls_map{lev} &  cls_map{lev} ~= ori));
            pyr_ttest_bc(ii, ori, lev) = p;
            [h p] = ttest2(temp(~cls_map{lev}), temp(cls_map{lev} == ori));
            pyr_ttest_ac(ii, ori, lev) = p;
            
        end
    end
    
    %Generate oriented CLS map
    [cls_map] =  mb_cls_map_orientations(cls, orientations);
    
    %Load pyramid
    dt = u_load(['C:\isbe\dev\background\dual_tree\mass_2\mass', zerostr(ii,3) ,'_dual_tree.mat']);
    
    %Compute values for each pixel typ in each sub-band
    for lev = 1:3
        
        cls_lev = imresize(cls_map{lev}, size(dt{lev}(:,:,1)), 'nearest');
        
        for ori = 1:6           
            
%             temp = ((dt{lev}(:,:,ori)));
%             %Compute means for each type
%             dt_all(ii, ori, lev) = mean(temp(:));
%             dt_all_cls(ii, ori, lev) = mean(temp(cls_lev>0));
%             dt_non_cls(ii, ori, lev) = mean(temp(~cls_lev));
%             dt_aligned_cls(ii, ori, lev) = mean(temp(cls_lev == ori+1));
%             dt_nonaligned_cls(ii, ori, lev) = mean(temp(cls_lev &  cls_lev ~= ori+1));
%             dt_all_nonaligned_cls(ii, ori, lev) = mean(temp(cls_lev ~= ori+1));
%             
%             temp = (abs(dt{lev}(:,:,ori)));
%             %Compute means for each type
%             dt_abs_all(ii, ori, lev) = mean(temp(:));
%             dt_abs_all_cls(ii, ori, lev) = mean(temp(cls_lev>0));
%             dt_abs_non_cls(ii, ori, lev) = mean(temp(~cls_lev));
%             dt_abs_aligned_cls(ii, ori, lev) = mean(temp(cls_lev == ori+1));
%             dt_abs_nonaligned_cls(ii, ori, lev) = mean(temp(cls_lev &  cls_lev ~= ori+1));
%             dt_abs_all_nonaligned_cls(ii, ori, lev) = mean(temp(cls_lev ~= ori+1));
            
            temp = (real(dt{lev}(:,:,ori)));
            [h p] = ttest2(temp(~cls_lev), temp(cls_lev &  cls_lev ~= ori+1));
            dt_ttest_real_ab(ii, ori, lev) = p;
            [h p] = ttest2(temp(cls_lev == ori+1), temp(cls_lev &  cls_lev ~= ori+1));
            dt_ttest_real_bc(ii, ori, lev) = p;
            [h p] = ttest2(temp(~cls_lev), temp(cls_lev == ori+1));
            dt_ttest_real_ac(ii, ori, lev) = p;
            
            temp = (imag(dt{lev}(:,:,ori)));
            [h p] = ttest2(temp(~cls_lev), temp(cls_lev &  cls_lev ~= ori+1));
            dt_ttest_imag_ab(ii, ori, lev) = p;
            [h p] = ttest2(temp(cls_lev == ori+1), temp(cls_lev &  cls_lev ~= ori+1));
            dt_ttest_imag_bc(ii, ori, lev) = p;
            [h p] = ttest2(temp(~cls_lev), temp(cls_lev == ori+1));
            dt_ttest_imag_ac(ii, ori, lev) = p;
            
        end
    end
    
    clear cls cls_map cls_lev dt pyr temp
end
%%
pyr_data = [...
    pyr_all(:)...
    pyr_all_cls(:)...
    pyr_non_cls(:)...
    pyr_aligned_cls(:)...
    pyr_nonaligned_cls(:)...
    pyr_all_nonaligned_cls(:) ];

dt_data = [...
    dt_all(:)...
    dt_all_cls(:)...
    dt_non_cls(:)...
    dt_aligned_cls(:)...
    dt_nonaligned_cls(:)...
    dt_all_nonaligned_cls(:) ];

dt_abs_data = [...
    dt_abs_all(:)...
    dt_abs_all_cls(:)...
    dt_abs_non_cls(:)...
    dt_abs_aligned_cls(:)...
    dt_abs_nonaligned_cls(:)...
    dt_abs_all_nonaligned_cls(:) ];
%%
pyr_data_sort = sortrows(pyr_data, 4);
pyr_data_sort_abs = sortrows(abs(pyr_data), 4);

dt_data_sort_abs = sortrows(abs(dt_data), 4);
dt_data_sort_real = sortrows(real(dt_data), 4);
dt_data_sort_imag = sortrows(imag(dt_data), 4);

dt_data_sort_real_abs = sortrows(abs(real(dt_data)), 4);
dt_data_sort_imag_abs = sortrows(abs(imag(dt_data)), 4);

dt_abs_data_sort = sortrows(dt_abs_data, 4);
%%
[dummy, pixel_class_pyr] = max(pyr_data_sort, [], 2);
[dummy, pixel_class_pyr_abs] = max(pyr_data_sort_abs, [], 2);
[dummy, pixel_class_dt_abs] = max(dt_data_sort_abs, [], 2);
[dummy, pixel_class_dt_real_abs] = max(dt_data_sort_real_abs, [], 2);
[dummy, pixel_class_dt_imag_abs] = max(dt_data_sort_imag_abs, [], 2);
[dummy, pixel_class_abs_dt] = max(dt_abs_data_sort, [], 2);

display(sum(pixel_class_pyr == 4) / length(pixel_class_pyr));
display(sum(pixel_class_pyr_abs == 4) / length(pixel_class_pyr_abs));
display(sum(pixel_class_dt_abs == 4) / length(pixel_class_dt_abs));
display(sum(pixel_class_dt_real_abs == 4) / length(pixel_class_dt_real_abs));
display(sum(pixel_class_dt_imag_abs == 4) / length(pixel_class_dt_imag_abs));
display(sum(pixel_class_abs_dt == 4) / length(pixel_class_dt_abs));
%%
figure; plot(pyr_data_sort, 'x');
figure; plot(pyr_data_sort_abs, 'x');
figure; plot(dt_data_sort_abs, 'x');
figure; plot(dt_data_sort_real, 'x');
figure; plot(dt_data_sort_imag, 'x');
figure; plot(dt_abs_data_sort, 'x');
%%
for xx = 0:10:40

orientations = [Inf -Inf; 0 pi; pi 2*pi; 2*pi 3*pi; -3*pi -2*pi; -2*pi -pi; -pi 0] / 6;

figure; hold on;
plot([0 41], [0 0], 'm:');
title(['Imaginary Coefficients Dual Tree - regions ', num2str(xx+1), ' to ', num2str(xx+10)]);

for ii = 1:10 %53 = number of masses
    %Load CLS result for each mass
    cls = u_load(['C:\isbe\dev\background\cls\mass_2\mass_cls_', zerostr(ii+xx,3), '.mat']);
    
    
    %Generate oriented CLS map
    [cls_map] =  mb_cls_map_orientations(cls, orientations);
    
    %Load pyramid
    dt = u_load(['C:\isbe\dev\background\dual_tree\mass_2\mass', zerostr(ii+xx,3) ,'_dual_tree.mat']);
    
    %Compute values for each pixel typ in each sub-band
    
    
    temp1 = [];
    temp2 = [];
    temp3 = [];
    
    for lev = 1:3
        
        cls_lev = imresize(cls_map{lev}, size(dt{lev}(:,:,1)), 'nearest');
        
        
        for ori = 1:6           
            
            temp = (imag(dt{lev}(:,:,ori)));
            
            
            temp1 = [temp1; temp(~cls_lev)];
            temp2 = [temp2; temp(cls_lev &  cls_lev ~= ori+1)];
            temp3 = [temp3; temp(cls_lev == ori+1)];
            
            plot((ii-1)*4 + 1, mean(temp(~cls_lev)), 'gx');
            plot((ii-1)*4 + 2, mean(temp(cls_lev &  cls_lev ~= ori+1)), 'gx');
            plot((ii-1)*4 + 3, mean(temp(cls_lev == ori+1)), 'gx');
        end
        
    end
    
    temp4 = NaN(max([length(temp1) length(temp2) length(temp3)]), 3);
    temp4(1:length(temp1), 1) = temp1;
    temp4(1:length(temp2), 2) = temp2;
    temp4(1:length(temp3), 3) = temp3;
    
    boxplot(temp4, 'positions', ((ii-1)*4 + 1):((ii-1)*4 + 3));
    
    clear cls cls_map cls_lev dt pyr temp*
end

axis([0 41 -10 10]);

end
%%
for xx = 0:10:40

figure; hold on;
plot([0 41], [0 0], 'm:');
title(['Pyramid Coefficients - regions ', num2str(xx+1), ' to ', num2str(xx+10)]);

for ii = 1:10 %53 = number of masses
    %Load CLS result for each mass
    cls = u_load(['C:\isbe\dev\background\cls\mass_2\mass_cls_', zerostr(ii+xx,3), '.mat']);
    
    
    %Generate oriented CLS map
    [cls_map] =  mb_cls_map_orientations(cls);
    
    %Load pyramid
    pyr = u_load(['C:\isbe\dev\background\pyramid5\mass_2\mass', zerostr(ii+xx,3) ,'_pyramid.mat']);
    
    %Compute values for each pixel typ in each sub-band
    
    
    temp1 = [];
    temp2 = [];
    temp3 = [];
    
    for lev = 1:3
        cls_lev = cls_map{lev};
        
        for ori = 1:5           
            
            temp = pyr{lev+1, ori};
            
            
            temp1 = [temp1; temp(~cls_lev)];
            temp2 = [temp2; temp(cls_lev &  cls_lev ~= ori+1)];
            temp3 = [temp3; temp(cls_lev == ori+1)];
            
            plot((ii-1)*4 + 1, mean(temp(~cls_lev)), 'gx');
            plot((ii-1)*4 + 2, mean(temp(cls_lev &  cls_lev ~= ori+1)), 'gx');
            plot((ii-1)*4 + 3, mean(temp(cls_lev == ori+1)), 'gx');
        end
        
    end
    
    temp4 = NaN(max([length(temp1) length(temp2) length(temp3)]), 3);
    temp4(1:length(temp1), 1) = temp1;
    temp4(1:length(temp2), 2) = temp2;
    temp4(1:length(temp3), 3) = temp3;
    
    boxplot(temp4, 'positions', ((ii-1)*4 + 1):((ii-1)*4 + 3));
    
    clear cls cls_map cls_lev dt pyr temp
end

axis([0 41 -10 10]);

end
%%    
for lev = 1:3
    for ori = 1:5           

        %Display box-plots for each type in each in each sub-band
        
        figure; boxplot(([...
            mean_all(:, ori, lev)...
            mean_all_cls(:, ori, lev)...
            mean_non_cls(:, ori, lev)...
            mean_aligned_cls(:, ori, lev)...
            mean_nonaligned_cls(:, ori, lev)...
            mean_all_nonaligned_cls(:, ori, lev) ]));

    end
end
%%
figure; boxplot(([...
    mean_all(:)...
    mean_all_cls(:)...
    mean_non_cls(:)...
    mean_aligned_cls(:)...
    mean_nonaligned_cls(:)...
    mean_all_nonaligned_cls(:) ]));
%%    
for lev = 1:3
    for ori = 1:6           

        %Display box-plots for each type in each in each sub-band
        
        figure; boxplot(([...
            mean_all_dt(:, ori, lev)...
            mean_all_cls_dt(:, ori, lev)...
            mean_non_cls_dt(:, ori, lev)...
            mean_aligned_cls_dt(:, ori, lev)...
            mean_nonaligned_cls_dt(:, ori, lev)...
            mean_all_nonaligned_cls_dt(:, ori, lev) ]));

    end
end
%%
figure; boxplot(([...
    mean_all_dt(:)...
    mean_all_cls_dt(:)...
    mean_non_cls_dt(:)...
    mean_aligned_cls_dt(:)...
    mean_nonaligned_cls_dt(:)...
    mean_all_nonaligned_cls_dt(:) ]));
%%
figure; boxplot(([...
    mean_all_dta(:)...
    mean_all_cls_dta(:)...
    mean_non_cls_dta(:)...
    mean_aligned_cls_dta(:)...
    mean_nonaligned_cls_dta(:)...
    mean_all_nonaligned_cls_dta(:) ]));
%%
pyr = u_load('C:\isbe\dev\background\pyramid5\mass_2\mass006_pyramid.mat');
figure; colormap(jet(256)); mask = ~roipoly(pyr{3,1});
%%
[synthesised_image_3_10 cluster_image_3_10] = mb_gmm_tex_synthesis(...
    'PathToTextureGMM', 'C:\isbe\dev\background\models\cls\mass_2_cls_k10_c_model_3_1_0',...
      'TargetImage', pyr{3,1},...
      'FilledImage', mask);

%%
cls = u_load('C:\isbe\dev\background\cls\mass_2\mass_cls_006.mat');
cls_map = mb_cls_map_orientations(cls);
cls_mask = cls_map{2} == 1;
%%
temp_band = pyr{3,1};
temp_band(~mask) = NaN;
%%
[synthesised_image_3_11 cluster_image_3_11] = mb_gmm_tex_synthesis_cls(...
    'PathToNormalGMM', 'C:\isbe\dev\background\models\cls\mass_2_cls_k10_c_model_3_1_0',...
    'PathToClsGMM', 'C:\isbe\dev\background\models\cls\mass_2_cls_k10_c_model_3_1_1',...
    'TargetImage', temp_band,...
    'FilledImage', mask | ~cls_mask,...
    'ClsMap', cls_mask,...
    'WindowSize2', 0,...
    'SampleImage2', [],...
    'ClusterImage', cluster_map,...
    'ForceCluster', 0,...
    'SynthesisMode', 'patch');
%%
filled_image = imread('C:\isbe\dev\background\misc\filled_image_3_1.bmp');
cls_map = imread('C:\isbe\dev\background\misc\cls_map_3_1.bmp');
%%
[synthesised_image_3_1 cluster_image_3_1] = mb_gmm_tex_synthesis_cls(...
    'PathToNormalGMM', 'C:\isbe\dev\background\models5\cls\mass_2_cls_k10_c_model_3_1_0',...
    'PathToClsGMM', 'C:\isbe\dev\background\models5\cls\mass_2_cls_k10_c_model_3_1_1',...
    'TargetImage', pyr{3,1},...
    'FilledImage', filled_image,...
    'ClsMap', cls_map,...
    'WindowSize2', 0,...
    'SampleImage2', [],...
    'SynthesisMode', 'pixel',...
    'MovieFilename', 'C:\isbe\epsrc_meeting\movies\temp.gif');
%%
[synthesised_image_3_2 cluster_image_3_2] = mb_gmm_tex_synthesis_cls(...
    'PathToNormalGMM', 'C:\isbe\dev\background\models5\cls\mass_2_cls_k10_c_model_3_2_0',...
    'PathToClsGMM', 'C:\isbe\dev\background\models5\cls\mass_2_cls_k10_c_model_3_2_1',...
    'TargetImage', pyr{3,2},...
    'FilledImage', filled_image,...
    'ClsMap', filled_image,...
    'WindowSize2', 0,...
    'SampleImage2', [],...
    'SynthesisMode', 'pixel');
[synthesised_image_3_3 cluster_image_3_3] = mb_gmm_tex_synthesis_cls(...
    'PathToNormalGMM', 'C:\isbe\dev\background\models5\cls\mass_2_cls_k10_c_model_3_3_0',...
    'PathToClsGMM', 'C:\isbe\dev\background\models5\cls\mass_2_cls_k10_c_model_3_3_1',...
    'TargetImage', pyr{3,3},...
    'FilledImage', filled_image,...
    'ClsMap', filled_image,...
    'WindowSize2', 0,...
    'SampleImage2', [],...
    'SynthesisMode', 'pixel');
[synthesised_image_3_4 cluster_image_3_4] = mb_gmm_tex_synthesis_cls(...
    'PathToNormalGMM', 'C:\isbe\dev\background\models5\cls\mass_2_cls_k10_c_model_3_4_0',...
    'PathToClsGMM', 'C:\isbe\dev\background\models5\cls\mass_2_cls_k10_c_model_3_4_1',...
    'TargetImage', pyr{3,4},...
    'FilledImage', filled_image,...
    'ClsMap', filled_image,...
    'WindowSize2', 0,...
    'SampleImage2', [],...
    'SynthesisMode', 'pixel');
[synthesised_image_3_5 cluster_image_3_5] = mb_gmm_tex_synthesis_cls(...
    'PathToNormalGMM', 'C:\isbe\dev\background\models5\normal_2\normal_2_a_model_3_1',...
    'PathToClsGMM', 'C:\isbe\dev\background\models5\cls\mass_2_cls_k10_c_model_3_1_1',...
    'TargetImage', pyr{3,5},...
    'FilledImage', filled_image,...
    'ClsMap', filled_image,...
    'WindowSize2', 0,...
    'SampleImage2', [],...
    'SynthesisMode', 'pixel');
%%
syn_args.FilledImage = filled_image;
syn_args.ModelDir = 'C:\isbe\dev\background\models5\cls\';
syn_args.ModelName = 'mass_2_cls_k10_c_model';
syn_args.ModelNameCls = 'mass_2_cls_k10_c_model';
syn_args.ClsMap = cls_map;
syn_args.TargetPyramid = pyr;
syn_args.CutOffLevel = 3;

[synthesised_image, pyramid, cluster_image] = mb_gmm_pyr_synthesis_cls(syn_args);
%%
%test function for randomly selecting regions
mb_sample_valid_region('MammoDir', 'C:\isbe\mammograms\new_CAD\BMP_2004\',...
    'MammoFiles', normal_list, ...
    'SaveDir', 'C:\isbe\mammograms\new_CAD\BMP_2004\temp\',...
    'WindowSize', 1024);


