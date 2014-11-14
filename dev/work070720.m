for ii = 1:20
    
    load(['C:\isbe\dev\annotations\', unchanged(ii).name]);
    figure('WindowStyle', 'docked');
    subplot(1,2,1)
    imagesc(double(mass.mass_ROI) - mass.mass_sub_smooth); axis image;
    
    load(['C:\isbe\dev\old_annotations\', unchanged(ii).name]);
    subplot(1,2,2)
    imagesc(double(mass.mass_ROI) - mass.mass_sub_smooth); axis image;
end