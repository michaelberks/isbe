mam_names = u_load('C:\isbe\asymmetry_project\data\mam_names\2004_screening_abnormals.mat');

roi_dir = 'C:\isbe\asymmetry_project\data\mammograms\2004_screening_processed\mass_roi\';
for ii = 1:10
    
    roi = u_load([roi_dir mam_names{ii} '_roi.mat']);
    
    line_map1 = load_uint8(['C:\isbe\asymmetry_project\data\line_maps\rf_prob\2004_screening_processed\mass_roi\'...
        mam_names{ii} '_roi.mat']);
    line_map2 = load_uint8([roi_dir 'results\13835\' mam_names{ii} '_roi_class.mat']);
    line_map3 = load_uint8([roi_dir 'results\13925\' mam_names{ii} '_roi_class.mat']);
    
    figure;
    subplot(2,2,1); imagesc(roi); axis image; colormap(gray(256));
    subplot(2,2,2); imagesc(line_map1); axis image; colormap(gray(256));
    subplot(2,2,3); imagesc(line_map2); axis image; colormap(gray(256));
    subplot(2,2,4); imagesc(line_map3); axis image; colormap(gray(256));
end
%%
ii = 1;
jj = 0;
while ii <= 50
    jj = jj + 1;
    if ismember(jj, [4 5 6 7 9 11]); continue; end
    
    copyfile(...
        ['C:\isbe\ug_project_options\2011\chitra\blobs1\blob' zerostr(jj,3) '_a.mat'],...
        ['C:\isbe\ug_project_options\2011\chitra\blobs\blob' zerostr(ii,3) '_a.mat']);
    copyfile(...
        ['C:\isbe\ug_project_options\2011\chitra\blobs1\blob' zerostr(jj,3) '_b.mat'],...
        ['C:\isbe\ug_project_options\2011\chitra\blobs\blob' zerostr(ii,3) '_b.mat']);
    
    copyfile(...
        ['C:\isbe\ug_project_options\2011\chitra\blobs1\blob' zerostr(jj,3) '_a.jpg'],...
        ['C:\isbe\ug_project_options\2011\chitra\blobs\blob' zerostr(ii,3) '_a.jpg']);
    copyfile(...
        ['C:\isbe\ug_project_options\2011\chitra\blobs1\blob' zerostr(jj,3) '_b.jpg'],...
        ['C:\isbe\ug_project_options\2011\chitra\blobs\blob' zerostr(ii,3) '_b.jpg']);
    
    ii = ii + 1;
end
     