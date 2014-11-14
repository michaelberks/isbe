function [roi] = segment_breast_mono(im_in, min_wav, mult, sigma);
%
% author: Michael Berks
% date: 04/08/2006 16:04
% function: segment breast using monogenic signal method, (Bo & Brady 2006)

if nargin == 1
	min_wav = 4;  % Default values for monogenic signal
    mult = 4;
    sigma = 0.41;
end
i1 = im_in(50:end-50, :);
%compute monogenic signal and local amplitude, phase and orientation
[local_amp, local_ori, local_phase] = monogenic(i1, 3, min_wav, mult, sigma, 1);
[rows cols] = size(i1);

%figure, imagesc(local_phase{2}); axis image;

phase_mask1 = local_phase{4} > pi/2;
phase_mask2 = local_phase{3} > pi/2;

ori_map(1:rows, 1:cols) = 1;
ori_map(local_ori{1} > -pi/2) = false;%1;
%ori_map(local_ori{1} > 0) = 2;
ori_map(local_ori{1} > 3*pi/4) = 1;%3;

%ori_smooth = medfilt2(ori_map, [5, 5]);
%ori_smooth = imopen(ori_map, strel('disk', 10, 0));
ori_smooth = imfill(ori_map);
[labeled, numObjects] = bwlabel(ori_smooth, 4);
ori_data = regionprops(labeled,'Area', 'PixelIdxList');
main_area = ori_data(find([ori_data.Area] == max([ori_data.Area]))).PixelIdxList;
ori_mask(rows, cols) = false;
ori_mask(main_area) = 1;

%morphological smoothing
%phase_mask = imopen(phase_mask, strel('disk', 10, 0));
%figure, imshow(phase_mask);
if 0
[start_r start_c] = find(phase_mask, 1);
%end_r = find(phase_mask(:, start_c), 1, 'last');

inner_border = bwtraceboundary(phase_mask, [start_r, start_c], 'E');

start_r = 1; start_c = find(ori_mask(1,:), 1, 'last');
outer_border = bwtraceboundary(ori_mask, [start_r, start_c], 'S');

end_r = rows; end_c = find(ori_mask(end_r, :), 1, 'last');
end_idx = find((outer_border(:,1) == end_r) & (outer_border(:,2) == end_c), 1);

outer_border = outer_border(1:end_idx, :);
end

roi = phase_mask1 & phase_mask2 & ori_mask;
[labeled, numObjects] = bwlabel(roi, 4);
roi_data = regionprops(labeled,'Area', 'PixelIdxList');
main_area = roi_data(find([roi_data.Area] == max([roi_data.Area]))).PixelIdxList;
roi_mask(rows, cols) = false;
roi_mask(main_area) = 1;

roi_phase1 = local_phase{4}(roi_mask);
%roi_phase2 = local_phase{3}(roi); 

roi_phase2 = local_phase{3}; roi_phase2(not(roi_mask)) = 0;
[max_vals max_idx] = max(roi_phase2');
breast_border = [1:rows; max_idx]';
figure, imagesc(roi_phase2); axis image;

if 0
temp1 = local_phase{3}; temp1(not(roi)) = 0;
temp2 = local_phase{4}; temp2(not(roi)) = 0;
figure, imagesc(temp1); axis image;
figure, imagesc(temp2); axis image;

roi_phase_round = (round(50*roi_phase));
phase_subs = (sort(unique(roi_phase_round)))/50;
aa = accumarray(roi_phase_round, roi_amp, [], @mean);
aa_idx = find(aa, 1);
%figure, plot(phase_subs, aa(aa_idx:end), 'x');
%hold on;

phase_idx = find(aa == max(aa)) - aa_idx + 1;
phase_thresh = phase_subs(phase_idx);
plot([phase_thresh phase_thresh], [0 max(aa)], 'r');

%breast_mask = local_phase{4} <= phase_thresh;
breast_mask = (local_phase{3} > pi-.1) & roi;
[breast_y breast_x] = find(breast_mask);
%hysteris thresholding of local_phase{2}
%aboveT2 = im > T2;                     % Edge points above lower
%                                           % threshold. 
%[aboveT1r, aboveT1c] = find(im > T1);  % Row and colum coords of points
                                       % above upper threshold.

% Obtain all connected regions in aboveT2 that include a point that has a
% value above T1 
%bw = bwselect(aboveT2, aboveT1c, aboveT1r, 8);


[start_r start_c] = find(breast_mask, 1);
breast_border = bwtraceboundary(breast_mask, [start_r, start_c], 'E');
end

figure, imagesc(im_in); colormap gray; axis image;
hold on;
%plot(breast_x, breast_y, 'r.');
plot(breast_border(:,2), 50+breast_border(:,1), 'r');

figure, imagesc(local_phase{3}); axis image;
hold on;
%plot(breast_x, breast_y, 'g.');
plot(breast_border(:,2), breast_border(:,1), 'g');

if 0;
%plot the initial breast border over the contrast enhanced image

figure, imagesc(im_in); colormap gray; axis image;
hold on;
plot(inner_border(:,2), inner_border(:,1), 'r');

figure, imagesc(ori_mask); axis image;
hold on;
plot(inner_border(:,2), inner_border(:,1), 'g');
plot(outer_border(:,2), outer_border(:,1), 'y');
end