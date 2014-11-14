 for ii = 1:2
    
    s1 = load(['M:\asymmetry_project\data\synthetic_lines\curves512\image' num2str(ii) '.mat']);
    s2 = load(['M:\asymmetry_project\data\synthetic_lines\match512\image' num2str(ii) '.mat']);
    
    figure;
    subplot(1,2,1); imagesc(s1.test_image); axis image;
    subplot(1,2,2); imagesc(s2.test_image); axis image;
end
%%
for ii = 1:2
    
    p1 = u_load(['Z:\asymmetry_project\data\synthetic_lines\circles512\results\191630\probability_image' zerostr(ii,3) '.mat']);
    p2 = u_load(['Z:\asymmetry_project\data\synthetic_lines\lines512\results\191630\probability_image' zerostr(ii,3) '.mat']);
    
    figure;
    subplot(1,2,1); imagesc(p1); axis image;
    subplot(1,2,2); imagesc(p2); axis image;
end
%%
copy_rf_from_hydra(191658);
%%
