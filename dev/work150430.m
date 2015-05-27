for i_xdim = 1:12

    [x_sort sort_idx] = sort(sample_data(:,i_xdim));

    ori_sort = ori_data(sort_idx);
    ori_c = complex(cos(ori_sort), sin(ori_sort));
    num_pts = length(ori_c);

    ssq = zeros(num_pts-1,1);
    dispersion = zeros(num_pts-1,1);
    ssq_lin = zeros(num_pts-1,1);
    for i_pt = 1:num_pts-1

        ori_l = ori_c(1:i_pt);
        ori_r = ori_c(i_pt+1:end);
        
        ang_l = ori_sort(1:i_pt);
        ang_r = ori_sort(i_pt+1:end);

        ssq(i_pt) = (abs(mean(ori_l)) + abs(mean(ori_r)));
        dispersion(i_pt) = ssq(i_pt) - abs(mean(ori_l) + mean(ori_r));
        ssq_lin(i_pt) = (sum(ang_l)^2)/i_pt + (sum(ang_r)^2)/(num_pts-i_pt);
    end
    
    [min_val, min_idx] = max(dispersion);
    ori_l = ori_c(1:min_idx);
    ori_r = ori_c(min_idx+1:end);

    ori_l = atan2(imag(ori_l), real(ori_l));
    ori_r = atan2(imag(ori_r), real(ori_r));
    
    counts_l = hist(ori_l, linspace(-pi,pi,24));
    counts_r = hist(ori_r, linspace(-pi,pi,24));
    
    figure;
    subplot(2,2,1); plot(ssq);
    subplot(2,2,2); plot(dispersion);
    subplot(2,2,3); plot(ssq_lin);
    subplot(2,2,4); hold on;
    bar(linspace(-pi,pi,24), counts_l, 0.5, 'b');
    bar(linspace(-pi,pi,24) + pi/24, counts_r, 0.5, 'r');
    title(['Min value = ' num2str(min_val)]);
    
end
%%
X = sample_data(:,13:24);
y = complex(sample_data(:,61), sample_data(:,62));

[tree assignednode] = tree_reg_train(X, y,...
             'random_m', 12,...
             'split_criterion', 'bob',...
			 'var_criterion', 'bob',...
             'split_min', 10, ...
             'end_cut_min', 0, ...
			 'w_prior', 0, ...
             'prune', 0,...
             'do_circular', [],...
             'do_ubound', 0,...
             'impure_thresh', 1e-4,...
             'names', []);
         
%%
pred_d = load('C:/isbe/nailfold/images/BII_comp/10598c_cxx_d_prediction.txt'); 
pred_o = read_complex_txt('C:/isbe/nailfold/images/BII_comp/10598c_cxx_o_prediction.txt'); 
pred_w = load('C:/isbe/nailfold/images/BII_comp/10598c_cxx_w_prediction.txt');         