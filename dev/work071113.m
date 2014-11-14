load C:\isbe\dev\files\u_files.mat

for ii = 21:53
    load(['C:\isbe\dev\masses\', bg_files(ii).name]);
    figure;
    subplot(1,2,1);
    imagesc(mass.mass_ROI); axis image; colormap(gray(256));
    subplot(1,2,2);
    imagesc(double(mass.mass_ROI) - double(mass.subtract_ROI));
    axis image; colormap(gray(256)); hold on;
    plot(mass.mass_outline(:,1), mass.mass_outline(:,2), 'y:');
    clear mass;
end

%%
for ii = 1:20
    load(['C:\isbe\dev\masses\', bg_files(ii).name]);
    figure;
    imagesc(double(mass.mass_ROI) - double(mass.subtract_ROI));
    axis image; colormap(gray(256)); hold on;
%     plot(mass.mass_outline(:,1), mass.mass_outline(:,2), 'y:');
    clear mass;
end