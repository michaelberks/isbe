p1 = u_load('Z:\asymmetry_project\data\mammograms\2004_screening\normals\results\lines_01\001LML_classification.mat');
p2 = u_load('Z:\asymmetry_project\data\mammograms\2004_screening\normals\results\lines_02\001LML_classification.mat');
mammo = u_load('Z:\asymmetry_project\data\mammograms\2004_screening\normals\001LML.mat');
mask = ~(~p1);

figure;
a1 = subplot(1,2,1); imagesc(1-p1); axis image; colormap(jet(256));
a2 = subplot(1,2,2); imagesc(1-p2); axis image; colormap(jet(256));
linkaxes([a1 a2]);

figure; imagesc(p1-p2); axis image; colormap(jet(256)); caxis([-1 1]); colorbar;

figure; hist(p1(mask)-p2(mask), linspace(-1, 1, 100)); colormap jet;

figure;
a1 = subplot(1,2,1); imagesc(mammo); axis image; colormap(jet(256));
a2 = subplot(1,2,2); imagesc(imfilter(1-p2, fspecial('gaussian', 40, 8))); axis image; colormap(jet(256));
linkaxes([a1 a2]);
%%
p1 = u_load('Z:\asymmetry_project\data\mammograms\2004_screening\normals\results\roi_01\001LML_classification.mat');
p2 = u_load('Z:\asymmetry_project\data\mammograms\2004_screening\normals\results\roi_02\001LML_classification.mat');
mammo = u_load('Z:\asymmetry_project\data\mammograms\2004_screening\normals\001LML.mat');
mask = ~(~p1);

figure;
a1 = subplot(1,2,1); imagesc(1-p1); axis image; colormap(jet(256));
a2 = subplot(1,2,2); imagesc(1-p2); axis image; colormap(jet(256));
linkaxes([a1 a2]);

figure; imagesc(p1-p2); axis image; colormap(jet(256)); caxis([-1 1]); colorbar;

figure; hist(p1(mask)-p2(mask), linspace(-1, 1, 100)); colormap jet;

figure;
a1 = subplot(1,2,1); imagesc(mammo); axis image; colormap(jet(256));
a2 = subplot(1,2,2); imagesc(imfilter(1-p2, fspecial('gaussian', 40, 8))); axis image; colormap(jet(256));
linkaxes([a1 a2]);
%%

mam_list = dir('C:\isbe\mammograms\new_CAD\BMP_2001\*.bmp');
for ii = 211:223%length(mam_list)
    
    mam = imread(['C:\isbe\mammograms\new_CAD\BMP_2001\' mam_list(ii).name]);
    
    %Check mammogram is portrait orientated
    [r c] = size(mam);
    if c > r
        mam = rot90(mam);
    end

    % resize to have 1024 rows
    mam = imresize(mam, [1024 NaN], 'bilinear');
    %right = ~isempty(strfind(mam_list(ii).name, 'R'));
    [r c] = size(mam);
    
    r1 = round(r / 3);
    r2 = 2*r1;
    c_half = round(c / 2);
    
    mam_mean = mean(mam(:));
    centre_mean = mean(mam(r1:r2,:));
    
    sum1 = sum(centre_mean(1:c_half) > mam_mean);
    sum2 = sum(centre_mean(c_half+1:end) > mam_mean);

    figure; imagesc(mam); axis image; hold on;
    
    %if sum1 > sum2 then the breast is on the left of the image
    if (sum1 > sum2)
        plot(10, r/2, 'r*', 'MarkerSize', 10);
    else
        plot(c-10, r/2, 'r*', 'MarkerSize', 10);
    end
    title(mam_list(ii).name);
    
    
%     figure; hold on;
%     plot(1:c, centre_mean);
%     plot([1 c], [mam_mean mam_mean], 'r');
%     if (sum1 > sum2)
%         plot(10, mam_mean, 'g*', 'MarkerSize', 10);
%     else
%         plot(c-10, mam_mean, 'g*', 'MarkerSize', 10);
%     end
%     
%     title(mam_list(ii).name);
    
end
%%
mkdir Z:\asymmetry_project\data\line_maps\2004_screening\abnormals\contralateral
mkdir Z:\asymmetry_project\data\orientation_maps\2004_screening\abnormals\contralateral

data_list = dir('Z:\asymmetry_project\data\line_maps\2004_screening\abnormals\*.mat');
for ii = 1:length(data_list)
    
    if ~exist(['Z:\asymmetry_project\data\mammograms\2004_screening\abnormals\meta\' data_list(ii).name(1:6) '_meta.mat'], 'file')
        display(['Moving ' data_list(ii).name(1:6)]);
        
        movefile(...
            ['Z:\asymmetry_project\data\line_maps\2004_screening\abnormals\' data_list(ii).name],...
            ['Z:\asymmetry_project\data\line_maps\2004_screening\abnormals\contralateral\' data_list(ii).name]);
        
        movefile(...
            ['Z:\asymmetry_project\data\orientation_maps\2004_screening\abnormals\' data_list(ii).name],...
            ['Z:\asymmetry_project\data\orientation_maps\2004_screening\abnormals\contralateral\' data_list(ii).name]);
    end
end
%%
c_list = dir('Z:\asymmetry_project\data\mammograms\2004_screening\normals\results\roi_01\*.mat');
for ii = 3:length(c_list)
    load(['Z:\asymmetry_project\data\mammograms\2004_screening\normals\results\roi_01\' c_list(ii).name]);
    
    probability_image = 1 - probability_image;
    save(['Z:\asymmetry_project\data\mammograms\2004_screening\normals\results\roi_01\' c_list(ii).name], 'probability_image');
end
    
    