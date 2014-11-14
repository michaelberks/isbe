syn_means = zeros(100,1);
syn_stds = zeros(100,1);
syn_local_stds = zeros(100,1);
real_means = zeros(100,1);
real_stds = zeros(100,1);
real_local_stds = zeros(100,1);

for ii = 1:100
    
    syn_bg = u_load(['C:\isbe\asymmetry_project\data\synthetic_backgrounds\new512\test\bg' zerostr(ii,5) '.mat']);
    
    [syn_local_stds(ii) syn_means(ii)] = image_stats(syn_bg, 32);
    syn_stds(ii) = std(syn_bg(:));
    
    real_bg = u_load(['C:\isbe\asymmetry_project\data\synthetic_backgrounds\real512\test\bg' zerostr(ii,5) '.mat']);
    
    [real_local_stds(ii) real_means(ii)] = image_stats(real_bg, 32);
    real_stds(ii) = std(real_bg(:));
end
%%
figure; hist([real_means syn_means]);
figure; hist([real_stds syn_stds]);
figure; hist([real_local_stds syn_local_stds]);
%%
figure; hold on;
for ii = 1:size(o_errors,1);
    plot(ii, o_errors(ii,1), 'bo', 'markersize', 6, 'markerfacecolor', 'b');
    plot(ii, o_errors(ii,2), 'rx', 'markersize', 6);
    if ii == 1
        legend({'2009', '2010'}, 'location', 'southeast');
    end
    plot([ii ii], o_errors(ii,:), 'k');
    plot(ii, o_errors(ii,1), 'bo', 'markersize', 6, 'markerfacecolor', 'b');
    plot(ii, o_errors(ii,2), 'rx', 'markersize', 6);
    
end
%%
syn_means = zeros(100,1);
syn_stds = zeros(100,1);
for ii = 1:100
    
    syn_bg = u_load(['C:\isbe\asymmetry_project\data\synthetic_lines\real512\image' zerostr(ii,3) '.mat']);
    
    syn_means(ii) = mean(syn_bg(:));
    syn_stds(ii) = std(syn_bg(:));

end
%%
mkdir C:\isbe\asymmetry_project\data\misc\ori_labels
for ii = 1:100
    s = load(['C:\isbe\asymmetry_project\data\synthetic_lines\real512\labels\label' zerostr(ii,3) '.mat']);
    ori_map = s.label_orientation;
    ori_map(isnan(ori_map)) = 0;
    save(['C:\isbe\asymmetry_project\data\misc\ori_labels\label' zerostr(ii,3) '.mat'], 'ori_map');
end