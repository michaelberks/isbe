%Evaluate initial detection methods
nf_files = dir('C:\isbe\nailfold\images\anonymous_oct\annotations_qt\*.txt');

test_dir = 'C:\isbe\nailfold\data\training\';
im_dir = [test_dir 'images\'];
votes_dir = [test_dir 'predictions\apex_votes\'];
apex_dir = [test_dir 'apexes\'];

all_candidate_labels_cu = [];
all_candidate_labels_ca = [];
all_candidate_labels_rf = [];

all_candidate_vals_cu = [];
all_candidate_vals_ca = [];
all_candidate_vals_rf = [];

num_gt = 0;

for nn = 7:12;%1:6%7:12
    %load in mask
    load(['C:\isbe\nailfold\images\anonymous_oct\masks\' nf_files(nn).name(1:end-11) '_mask.mat'], 'nailfold_mask');
    
%     %load in the correlation masks - unaligned
%     corr_unaligned = load([nailfoldroot 'data\apex_detection\unaligned\' nf_files(nn).name(1:end-4) '_corr.mat'], 'C1', 'C2');
%     corr_unaligned = corr_unaligned.C1 .* corr_unaligned.C2;
%     corr_unaligned(~nailfold_mask) = 0;
%     [maxima_pos_cu, maxima_vals_cu] = local_image_maxima(corr_unaligned, 20, [], 0);    
%     clear corr_unaligned;
%     
    %load in the correlation masks - aligned
    corr_aligned = load([nailfoldroot 'data\apex_detection\aligned\' nf_files(nn).name(1:end-4) '_corr.mat'], 'C1', 'C2');
    corr_aligned = corr_aligned.C1 .* corr_aligned.C2;
    corr_aligned(~nailfold_mask) = 0;
    [maxima_pos_ca, maxima_vals_ca] = local_image_maxima(corr_aligned, 20, [], 0);
    clear corr_aligned;
    
%     load in the RF vote maps
%     vote_rf = u_load([votes_dir 'nailfold' zerostr(nn, 3) '_votes.mat']);  
%     vote_rf = imfilter(vote_rf, fspecial('gaussian', [10 10], 2));
%     [maxima_pos_rf, maxima_vals_rf] = local_image_maxima(vote_rf, 20, [], 0);
%     figure; imgray(vote_rf);
%     clear vote_rf;
    
    %Load in ground truth apex locations
    apexes = u_load(['C:\isbe\nailfold\data\apexes\' nf_files(nn).name(1:end-4) '_apexes.mat']);
    apex_gt = zeros(length(apexes),2);
    for ii = 1:length(apexes)
        apex_gt(ii,:) = apexes{ii}(16,:);
    end
    num_gt = num_gt + length(apexes);
    
%     %Workout whether candidate are hits or not
%     [candidate_labels_cu, gt_hits_cu] = evaluate_apex_candidates(apex_gt, maxima_pos_cu, 15);
%     all_candidate_labels_cu = [all_candidate_labels_cu; candidate_labels_cu]; %#ok
%     all_candidate_vals_cu = [all_candidate_vals_cu; maxima_vals_cu]; %#ok
%     save([nailfoldroot 'data\apex_detection\embs\test\nailfold' zerostr(nn,3) '_candidates_cu.mat'], 'maxima_*_cu', 'candidate_labels_cu');
%     
    [candidate_labels_ca] = evaluate_apex_candidates(apex_gt, maxima_pos_ca, 15);
    all_candidate_labels_ca = [all_candidate_labels_ca; candidate_labels_ca]; %#ok
    all_candidate_vals_ca = [all_candidate_vals_ca; maxima_vals_ca]; %#ok
    save([nailfoldroot 'data\apex_detection\embs\nailfold' zerostr(nn,3) '_candidates_ca.mat'], 'maxima_*_ca', 'candidate_labels_ca');
        
    
%     [candidate_labels_rf] = evaluate_apex_candidates(apex_gt, maxima_pos_rf, 15);    
%     all_candidate_labels_rf = [all_candidate_labels_rf; candidate_labels_rf]; %#ok
%     all_candidate_vals_rf = [all_candidate_vals_rf; maxima_vals_rf]; %#ok
%     save([nailfoldroot 'data\apex_detection\embs\nailfold' zerostr(nn,3) '_candidates_rf.mat'], 'maxima_*_rf', 'candidate_labels_rf');
        
end
%%
%save([nailfoldroot 'data\apex_detection\embs\all_candidates.mat'], 'all_candidate_*');

[~, idx_cu] = sort(all_candidate_vals_cu, 'descend');
tp_count_cu = cumsum(all_candidate_labels_cu(idx_cu));
fp_count_cu = (1:length(all_candidate_labels_cu))' - tp_count_cu;
sensitivity_cu = tp_count_cu / num_gt;
thresh_cu = all_candidate_vals_cu(idx_cu( sum(sensitivity_cu < 0.98)+1 ));
save([nailfoldroot 'data\apex_detection\embs\thresh_cu.mat'], 'thresh_cu');

% [~, idx_ca] = sort(all_candidate_vals_ca, 'descend');
% tp_count_ca = cumsum(all_candidate_labels_ca(idx_ca));
% fp_count_ca = (1:length(all_candidate_labels_ca))' - tp_count_ca;
% sensitivity_ca = tp_count_ca / num_gt;
% thresh_ca = all_candidate_vals_ca(idx_ca( sum(sensitivity_ca < 0.98)+1 ));
% save([nailfoldroot 'data\apex_detection\embs\thresh_ca.mat'], 'thresh_ca');

[~, idx_rf] = sort(all_candidate_vals_rf, 'descend');
tp_count_rf = cumsum(all_candidate_labels_rf(idx_rf));
fp_count_rf = (1:length(all_candidate_labels_rf))' - tp_count_rf;
sensitivity_rf = tp_count_rf / num_gt;
thresh_rf = all_candidate_vals_rf(idx_rf( sum(sensitivity_rf < 0.98) + 1 ));
save([nailfoldroot 'data\apex_detection\embs\thresh_rf.mat'], 'thresh_rf');

display(['Number or candidates using unaligned correlation templates: ' num2str(sum(all_candidate_vals_cu > thresh_cu))]);
% display(['Number or candidates using aligned correlation templates: ' num2str(sum(all_candidate_vals_ca > thresh_ca))]);
% display(['Number or candidates using RF regression: ' num2str(sum(all_candidate_vals_rf > thresh_rf))]);

% figure; hold all; title('Template matching - unaligned shapes');
% plot(log10(fp_count_cu), sensitivity_cu);
% plot(log10(fp_count_ca), sensitivity_ca);
% plot(log10(fp_count_rf), sensitivity_rf);
% 
% legend({'Corr Unaligned', 'Corr Aligned', 'RF'}, 'location', 'southeast');
% set(gca, 'ylim', [0 1]);
%%
%--------------------------------------------------------------------------
%% Generate data for classification task
nf_files = dir('C:\isbe\nailfold\images\anonymous_oct\annotations_qt\*.txt');

%% Unaligned correlation
mkdir([nailfoldroot 'data\apex_detection\embs\unaligned\candidate_patches\']);
mkdir([nailfoldroot 'data\apex_detection\embs\unaligned\candidate_data\']);

load([nailfoldroot 'data\apex_detection\embs\thresh_cu.mat'], 'thresh_cu');
apex_num = 1;
for nn = 1:6%7:12
    
    %Load candidate data
    load([nailfoldroot 'data\apex_detection\embs\nailfold' zerostr(nn,3) '_candidates_cu.mat']);
    
    %Discard all candidates below the threshold (computed to save 98% of
    %TPs)
    discard_idx = maxima_vals_cu < thresh_cu;
    maxima_pos_cu(discard_idx,:) = []; 
    maxima_vals_cu(discard_idx,:) = [];
    candidate_labels_cu(discard_idx,:) = [];
    num_candidates = size(maxima_pos_cu,1);
    
    %load in nailfold
    image_path = ['C:\isbe\nailfold\images\anonymous_oct\bmp\' nf_files(nn).name(1:end-11) '.bmp'];
    nailfold = imread(image_path);   
    nailfold = nailfold(:,:,1);
    [rows cols] = size(nailfold);
    
    %loop through remaining candidates, saving patch and and text file
    %giving the potential apex
    for i_c = 1:num_candidates
        
        sr = max(1, floor(maxima_pos_cu(i_c,2) - 100));
        er = min(rows, ceil(maxima_pos_cu(i_c,2) + 100));
        sc = max(1, floor(maxima_pos_cu(i_c,1) - 100));
        ec = min(cols, floor(maxima_pos_cu(i_c,1) + 100));

        %Sample patch from image
        image_patch = nailfold(sr:er, sc:ec);
        imwrite(image_patch, [nailfoldroot 'data\apex_detection\embs\unaligned\candidate_patches\apex' zerostr(apex_num, 4) '.png']);

        %Write out a pts file we can read in to VXL
        fid1 = fopen([nailfoldroot 'data\apex_detection\embs\unaligned\candidate_data\apex' zerostr(apex_num, 4) '.pts'], 'wt');
        fprintf(fid1, '%s \n', '{'); 
        fprintf(fid1,'%.2f %.2f %.2f \n', maxima_pos_cu(i_c,1) - sc, maxima_pos_cu(i_c,2) - sr, candidate_labels_cu(i_c));
        fprintf(fid1, '%s \n', '}');
        fprintf(fid1, 'nailfold: %s \n', image_path);
        fprintf(fid1, '%s %d \n', 'start_row:', sr);
        fprintf(fid1, '%s %d \n', 'start_col: ', sc);
        fclose(fid1);
        
        %Icrement the vessel count
        apex_num = apex_num + 1;
    end
        
end
%% Aligned correlation
mkdir([nailfoldroot 'data\apex_detection\embs\test\aligned\candidate_patches\']);
mkdir([nailfoldroot 'data\apex_detection\embs\test\aligned\candidate_data\']);

load([nailfoldroot 'data\apex_detection\embs\training\thresh_ca.mat'], 'thresh_ca');
apex_num = 1;
for nn = 7:12%1:6%
    
    %Load candidate data
    load([nailfoldroot 'data\apex_detection\embs\test\nailfold' zerostr(nn,3) '_candidates_ca.mat']);
    
    %Discard all candidates below the threshold (computed to save 98% of
    %TPs)
    discard_idx = maxima_vals_ca < thresh_ca;
    maxima_pos_ca(discard_idx,:) = []; 
    maxima_vals_ca(discard_idx,:) = [];
    candidate_labels_ca(discard_idx,:) = [];
    num_candidates = size(maxima_pos_ca,1);
    
    %load in nailfold
    image_path = ['C:\isbe\nailfold\images\anonymous_oct\bmp\' nf_files(nn).name(1:end-11) '.bmp'];
    nailfold = imread(image_path);   
    nailfold = nailfold(:,:,1);
    [rows cols] = size(nailfold);
    
    %loop through remaining candidates, saving patch and and text file
    %giving the potential apex
    for i_c = 1:num_candidates
        
        sr = max(1, floor(maxima_pos_ca(i_c,2) - 100));
        er = min(rows, ceil(maxima_pos_ca(i_c,2) + 100));
        sc = max(1, floor(maxima_pos_ca(i_c,1) - 100));
        ec = min(cols, floor(maxima_pos_ca(i_c,1) + 100));

        %Sample patch from image
        image_patch = nailfold(sr:er, sc:ec);
        imwrite(image_patch, [nailfoldroot 'data\apex_detection\embs\test\aligned\candidate_patches\apex' zerostr(apex_num, 4) '.png']);

        %Write out a pts file we can read in to VXL
        fid1 = fopen([nailfoldroot 'data\apex_detection\embs\test\aligned\candidate_data\apex' zerostr(apex_num, 4) '.pts'], 'wt');
        fprintf(fid1, '%s \n', '{'); 
        fprintf(fid1,'%.2f %.2f %.2f \n', maxima_pos_ca(i_c,1) - sc, maxima_pos_ca(i_c,2) - sr, candidate_labels_ca(i_c));
        fprintf(fid1, '%s \n', '}');
        fprintf(fid1, 'nailfold: %s \n', image_path);
        fprintf(fid1, '%s %d \n', 'start_row:', sr);
        fprintf(fid1, '%s %d \n', 'start_col: ', sc);
        fclose(fid1);
        
        %Icrement the vessel count
        apex_num = apex_num + 1;
    end
        
end

%% RF votes
mkdir([nailfoldroot 'data\apex_detection\embs\rf\candidate_patches\']);
mkdir([nailfoldroot 'data\apex_detection\embs\rf\candidate_data\']);

load([nailfoldroot 'data\apex_detection\embs\thresh_rf.mat'], 'thresh_rf');
apex_num = 1;
for nn = 1:6%7:12
    
    %Load candidate data
    load([nailfoldroot 'data\apex_detection\embs\nailfold' zerostr(nn,3) '_candidates_rf.mat']);
    
    %Discard all candidates below the threshold (computed to save 98% of
    %TPs)
    discard_idx = maxima_vals_rf < thresh_rf;
    maxima_pos_rf(discard_idx,:) = []; 
    maxima_vals_rf(discard_idx,:) = [];
    candidate_labels_rf(discard_idx,:) = [];
    num_candidates = size(maxima_pos_rf,1);
    
    %load in nailfold
    image_path = ['C:\isbe\nailfold\images\anonymous_oct\bmp\' nf_files(nn).name(1:end-11) '.bmp'];
    nailfold = imread(image_path);   
    nailfold = nailfold(:,:,1);
    [rows cols] = size(nailfold);
    
    %loop through remaining candidates, saving patch and and text file
    %giving the potential apex
    for i_c = 1:num_candidates
        
        sr = max(1, floor(maxima_pos_rf(i_c,2) - 100));
        er = min(rows, ceil(maxima_pos_rf(i_c,2) + 100));
        sc = max(1, floor(maxima_pos_rf(i_c,1) - 100));
        ec = min(cols, floor(maxima_pos_rf(i_c,1) + 100));

        %Sample patch from image
        image_patch = nailfold(sr:er, sc:ec);
        imwrite(image_patch, [nailfoldroot 'data\apex_detection\embs\rf\candidate_patches\apex' zerostr(apex_num, 4) '.png']);

        %Write out a pts file we can read in to VXL
        fid1 = fopen([nailfoldroot 'data\apex_detection\embs\rf\candidate_data\apex' zerostr(apex_num, 4) '.pts'], 'wt');
        fprintf(fid1, '%s \n', '{'); 
        fprintf(fid1,'%.2f %.2f %.2f \n', maxima_pos_rf(i_c,1) - sc, maxima_pos_rf(i_c,2) - sr, candidate_labels_rf(i_c));
        fprintf(fid1, '%s \n', '}');
        fprintf(fid1, 'nailfold: %s \n', image_path);
        fprintf(fid1, '%s %d \n', 'start_row:', sr);
        fprintf(fid1, '%s %d \n', 'start_col: ', sc);
        fclose(fid1);
        
        %Icrement the vessel count
        apex_num = apex_num + 1;
    end
        
end
%%
load([nailfoldroot 'data\apex_detection\embs\training\thresh_ca.mat']);
keep_idx = all_candidate_vals_ca > thresh_ca;

thresh_cand_vals_ca = all_candidate_vals_ca(keep_idx,:);
thresh_cand_labels_ca = all_candidate_labels_ca(keep_idx,:);

class_ca = load('C:\isbe\nailfold\data\apex_detection\embs\test\aligned\class_output.txt');

[roc_pts_rf auc_rf] = calculate_roc_curve(class_ca, thresh_cand_labels_ca);
[roc_pts_ca auc_ca] = calculate_roc_curve(thresh_cand_vals_ca, thresh_cand_labels_ca);

figure; axis equal; axis([0 1 0 1]); hold on;
plot(roc_pts_rf(:,1), roc_pts_rf(:,2), 'r');
plot(roc_pts_ca(:,1), roc_pts_ca(:,2), 'b');
%%

%%
for apex_num = 401:419
    filename = [nailfoldroot 'data\apex_detection\embs\training\aligned\candidate_data\apex' zerostr(apex_num, 4) '.pts'];
    fid = fopen(filename, 'r');
    frewind(fid);
    s = textscan(fid, '%s', 'commentstyle', '//'); s = s{1};
    fclose(fid);
    if ~str2num(s{4})
        apex = imread(...
            [nailfoldroot 'data\apex_detection\embs\training\aligned\candidate_patches\apex' zerostr(apex_num, 4) '.png']);
        x = str2num(s{2});
        y = str2num(s{3});
        figure; imgray(apex); plot(x, y, 'rx');
        title(['Apex ' zerostr(apex_num, 4)]);
    end
end
%%
swap_idx = [7 8 12 12 16 17 21 22 23 24 28 43 46 66 88 108 109 114 124 130 134 135 137 138 139 140 142 147 ...
    149 150 155 156 157 159 162 165 166 168 171 172 173 174 175 176 177 178 186 187 193 198 ...
    202 203 205 231 237 238 254 301 315 327 344 347 365 372 377 386 387 388 389 393 399 402 ...
    404 409 410];
%%
%% Aligned correlation
mkdir([nailfoldroot 'data\apex_detection\embs\training\aligned\candidate_patches\']);
mkdir([nailfoldroot 'data\apex_detection\embs\training\aligned\candidate_data\']);

load([nailfoldroot 'data\apex_detection\embs\training\thresh_ca.mat'], 'thresh_ca');
apex_num = 1;
for nn = 1:6%
    
    %Load candidate data
    load([nailfoldroot 'data\apex_detection\embs\training\nailfold' zerostr(nn,3) '_candidates_ca.mat']);
    
    %Discard all candidates below the threshold (computed to save 98% of
    %TPs)
    discard_idx = maxima_vals_ca < thresh_ca;
    maxima_pos_ca(discard_idx,:) = []; 
    maxima_vals_ca(discard_idx,:) = [];
    candidate_labels_ca(discard_idx,:) = [];
    num_candidates = size(maxima_pos_ca,1);
    
    %load in nailfold
    image_path = ['C:\isbe\nailfold\images\anonymous_oct\bmp\' nf_files(nn).name(1:end-11) '.bmp'];
    nailfold = imread(image_path);   
    nailfold = nailfold(:,:,1);
    [rows cols] = size(nailfold);
    
    %loop through remaining candidates, saving patch and and text file
    %giving the potential apex
    for i_c = 1:num_candidates
        
        sr = max(1, floor(maxima_pos_ca(i_c,2) - 100));
        er = min(rows, ceil(maxima_pos_ca(i_c,2) + 100));
        sc = max(1, floor(maxima_pos_ca(i_c,1) - 100));
        ec = min(cols, floor(maxima_pos_ca(i_c,1) + 100));

        %Sample patch from image
        image_patch = nailfold(sr:er, sc:ec);
        imwrite(image_patch, [nailfoldroot 'data\apex_detection\embs\training\aligned\candidate_patches\apex' zerostr(apex_num, 4) '.png']);

        %Write out a pts file we can read in to VXL
        fid1 = fopen([nailfoldroot 'data\apex_detection\embs\training\aligned\candidate_data\apex' zerostr(apex_num, 4) '.pts'], 'wt');
        fprintf(fid1, '%s \n', '{'); 
        if ismember(apex_num, swap_idx)
            fprintf(fid1,'%.2f %.2f %.2f \n', maxima_pos_ca(i_c,1) - sc, maxima_pos_ca(i_c,2) - sr, 1);
        else
            fprintf(fid1,'%.2f %.2f %.2f \n', maxima_pos_ca(i_c,1) - sc, maxima_pos_ca(i_c,2) - sr, candidate_labels_ca(i_c));
        end
        fprintf(fid1, '%s \n', '}');
        fprintf(fid1, 'nailfold: %s \n', image_path);
        fprintf(fid1, '%s %d \n', 'start_row:', sr);
        fprintf(fid1, '%s %d \n', 'start_col: ', sc);
        fclose(fid1);
        
        %Icrement the vessel count
        apex_num = apex_num + 1;
    end
        
end
   