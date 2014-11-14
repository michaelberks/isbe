%%
for ii = 61:80
    %if inc(ii) > 0.5
        load(['annotations\', files(ii).name]);
        figure('WindowStyle', 'docked');
        subplot(1,2,1);
        image(mass.mass_ROI); colormap(gray(256));
        axis('image');
        subplot(1,2,2);
        imagesc(uint8(mass.mass_sub_smooth)); colormap(gray(256));
        %image(double(mass.mass_ROI) - mass.mass_sub_smooth); colormap(gray(256));
        %hold on;
        %plot(mass.mass_outline([1:end,1],1), mass.mass_outline([1:end,1],2));
        axis('image'); %title(['landmarks diff = ', num2str(inc(ii))]);
        clear mass;
    %end
end
%%
bg_var(185, 1) = 0;
mass_var(185, 1) = 0;
for ii = 1:185
    load(['annotations\', files(ii).name]);
    [im_var im_mean var_map1] = image_stats(mass.mass_ROI, 11);
    bg_var(ii) = mean(var_map1(mass.mass_list20_40));
    mass_var(ii) = mean(var_map1(mass.mass_list));
    %[im_var im_mean var_map2] = image_stats(double(mass.mass_ROI) -
    %mass.mass_sub_smooth, 11);
end
%%
for ii = 1:length(files)
    load(['annotations\', files(ii).name]);
    mask = roipoly(mass.mass_ROI, mass.mass_outline(:,1), mass.mass_outline(:,2));
    for jj = 1:20
        mask = imdilate(mask, strel('disk', 1));
    end
    mask2 = mask;
    for jj = 1:40
        mask2 = imdilate(mask2, strel('disk', 1));
    end
    mass.mass_list20_40 = find(mask2 - mask);
    save(['annotations\', files(ii).name], 'mass');
    clear mass mask*
end
%%
for ii = [61 69 104 120 124 126 135 141 154 163]
    load(['annotations\', files(ii).name]);
    figure('WindowStyle', 'docked');
    image(mass.mass_ROI); colormap(gray(256));
    axis('image');
    clear mass;
end