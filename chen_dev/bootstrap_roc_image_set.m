function auc_cell = bootstrap_roc_image_set(image_dir, probability_dir, bootstrap_number)

[roc_pts, auc, tp_counts, fp_counts, t_counts, f_counts] = compute_roc_image_set(image_dir, probability_dir);

rand('twister', sum(100*clock));

% bootstrap_number = [20, 40, 80, 160, 320, 640];

auc_cell=cell(0);
% samples = 320;
for m=1:length(bootstrap_number);
    total_auc = [];
    for k=1:bootstrap_number(m)
        bootsrap_idx=ceil(100*rand(100,1));
        %Compute ROC points for complete set of data
        auc = compute_auc(bootsrap_idx, tp_counts, fp_counts, t_counts, f_counts);
        total_auc(k)=auc;
        
        % disp(k);
    end
auc_cell{m}=total_auc; 

end
return

function auc = compute_auc(bootsrap_idx, tp_counts, fp_counts, t_counts, f_counts)
tp_counts = tp_counts(bootsrap_idx, :);
fp_counts = fp_counts(bootsrap_idx, :);
t_counts = t_counts(bootsrap_idx);
f_counts = f_counts(bootsrap_idx);
roc_pts = [sum(fp_counts)' / sum(f_counts) sum(tp_counts)' / sum(t_counts)];

%Compute AUC for ROC curve
auc = sum( (roc_pts(2:end,1)-roc_pts(1:end-1,1)) .* roc_pts(1:end-1,2)) + ...
    0.5*sum( (roc_pts(2:end,1)-roc_pts(1:end-1,1)) .* (roc_pts(2:end,2)-roc_pts(1:end-1,2)) );
return
