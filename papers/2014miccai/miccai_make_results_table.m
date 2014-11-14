%%
nb = 1000;
num_images = size(co_occurrence_t,4);

pr_ab = zeros(4,2,nb);
re_ab = zeros(4,2,nb);
ac_ab = zeros(4,2,nb);
ka_ab = zeros(4,2,nb);
fa_ab = zeros(4,2,nb);

for i_b = 1:nb
    b_idx = ceil(num_images*rand(num_images,1));
    co = co_occurrence_t(:,:,:,b_idx);
    
    for ii = 1:4
        %
        switch ii;
            case 1;
                a_v_b = sum(sum(co,4),3);
            case 2;
                a_v_b = squeeze(sum(sum(co,4),1))';
            case 3;
                a_v_b = squeeze(sum(sum(co,4),2))';
            case 4;
                a_b_c = sum(co,4);
                a_b_c(2,:,:) = a_b_c(2,:,:) + a_b_c(3,:,:); a_b_c(3,:,:) = [];
                a_b_c(:,2,:) = a_b_c(:,2,:) + a_b_c(:,3,:); a_b_c(:,3,:) = [];
                a_b_c(:,:,2) = a_b_c(:,:,2) + a_b_c(:,:,3); a_b_c(:,:,3) = [];

                a_v_b = zeros(2,2);
                a_v_b(1,1) = a_b_c(1,1,1);
                a_v_b(1,2) = a_b_c(2,2,1);
                a_v_b(2,1) = a_b_c(1,1,2);
                a_v_b(2,2) = a_b_c(2,2,2);
        end
        if ii < 4
            a_v_b(:,2) = a_v_b(:,2) + a_v_b(:,3); a_v_b(:,3) = [];
            a_v_b(2,:) = a_v_b(2,:) + a_v_b(3,:); a_v_b(3,:) = [];
        end
        total_a = sum(a_v_b,1);
        total_b = sum(a_v_b,2);

        pr_ab(ii,1, i_b) = 100*a_v_b(2,2) / total_b(2);
        re_ab(ii,1, i_b) = 100*a_v_b(2,2) / total_a(2);
        ac_ab(ii,1, i_b) = 100*sum(diag(a_v_b))/sum(a_v_b(:));
        a_v_b(1,1) = 1e5;
        ka_ab(ii,1, i_b) = compute_kappa_statistic(a_v_b);
        
        fa_ab(ii,1, i_b) = (2*pr_ab(ii,1, i_b)*re_ab(ii,1, i_b)) / (pr_ab(ii,1, i_b) + re_ab(ii,1, i_b));
    end
    
    co = sum(co,4);
    co(1,:,:) = [];
    co(:,1,:) = [];
    co(:,:,1) = [];
        
    for ii = 1:4

        switch ii;
            case 1;
                a_v_b = sum(co,3);
            case 2;
                a_v_b = squeeze(sum(co,1))';
            case 3;
                a_v_b = squeeze(sum(co,2))';
            case 4;
                a_v_b = zeros(2,2);
                a_v_b(1,1) = co(1,1,1);
                a_v_b(1,2) = co(2,2,1);
                a_v_b(2,1) = co(1,1,2);
                a_v_b(2,2) = co(2,2,2);
        end

        total_a = sum(a_v_b,1);
        total_b = sum(a_v_b,2);

        pr_ab(ii,2, i_b) = 100*a_v_b(2,2) / total_b(2);
        re_ab(ii,2, i_b) = 100*a_v_b(2,2) / total_a(2);
        ac_ab(ii,2, i_b) = 100*sum(diag(a_v_b))/sum(a_v_b(:));
        ka_ab(ii,2, i_b) = compute_kappa_statistic(a_v_b);
        fa_ab(ii,2, i_b) = (2*pr_ab(ii,2, i_b)*re_ab(ii,2, i_b)) / (pr_ab(ii,2, i_b) + re_ab(ii,2, i_b));
    end
end

pr_ab_m = (mean(pr_ab,3));
pr_ab_s = (std(pr_ab,1,3));

re_ab_m = (mean(re_ab,3));
re_ab_s = (std(re_ab,1,3));

ac_ab_m = (mean(ac_ab,3));
ac_ab_s = (std(ac_ab,1,3));

ka_ab_m = (mean(ka_ab,3));
ka_ab_s = (std(ka_ab,1,3));

fa_ab_m = (mean(fa_ab,3));
fa_ab_s = (std(fa_ab,1,3));
%%

%Write a table from this
o_txt = {'$O_2$ v $O_1$', '$O_3$ v $O_1$', '$O_3$ v $O_2$', '$O_3$ v $O_1,O_2$'};

fid = fopen('C:\isbe\matlab_code\mab\papers\2014miccai\results_table.txt', 'wt');
fprintf(fid, '%s \n', '\begin{tabular*}{0.95\textwidth}{@{\extracolsep{\fill} } l r r r r r}');
fprintf(fid, '%s \n', '\toprule');
fprintf(fid, '%s \n', '%');

fprintf(fid, '%s \n', 'Observers		& \multicolumn{3}{c}{Capillary detection} 	& \multicolumn{2}{c}{Distal v Non-distal}	\\');
fprintf(fid, '%s \n', '			& Precision    		& Recall		& $F$-measure		& Accuracy		& Cohen''s $\kappa$ 	\\');
for ii = 1:4
    fprintf(fid, '%s', o_txt{ii});
    fprintf(fid, '& %3.1f $\\pm$ %2.1f ', pr_ab_m(ii,1), pr_ab_s(ii,1));
    fprintf(fid, '& %3.1f $\\pm$ %2.1f ', re_ab_m(ii,1), re_ab_s(ii,1));
    fprintf(fid, '& %3.1f $\\pm$ %2.1f ', fa_ab_m(ii,1), fa_ab_s(ii,1));
    fprintf(fid, '& %3.1f $\\pm$ %2.1f ', ac_ab_m(ii,2), ac_ab_s(ii,2));
    fprintf(fid, '& %4.3f $\\pm$ %3.3f ', ka_ab_m(ii,2), ka_ab_s(ii,2));
    fprintf(fid, '%s \n', '\\');
end

fprintf(fid, '%s \n', '%');
fprintf(fid, '%s \n', '\bottomrule \noalign{\smallskip}');
fprintf(fid, '%s \n', '\end{tabular*}');
fclose(fid);