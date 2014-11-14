%%
for ii = 101:185
    load(['annotations\', files(ii).name]);
    
    figure('WindowStyle', 'docked', 'Name', mass.name);
    im_large = imread(mass.name);
    %image(im_large); colormap(gray(256));
    %hold('on'); axis('image');
    outline_x = mass.mass_outline(:,1) + mass.C1 - 1;
    outline_y = mass.mass_outline(:,2) + mass.R1 - 1;
    %plot(outline_x([1:end, 1]), outline_y([1:end, 1]));
    
    %figure('WindowStyle', 'docked', 'Name', mass.name);
    im_small = imresize(im_large, 0.25, 'nearest'); clear im_large;
    im_smooth = imfilter(im_small, fspecial('gaussian', 15, 3), 'symmetric');
    image(im_smooth); colormap(gray(256));
    hold('on'); axis('image');
    plot(outline_x([1:end, 1])/4, outline_y([1:end, 1])/4);
    
    clear mass;
end
%%
for ii = [3 14 17 39 82 113 156 160 184]
    
    load(['..\..\dev\annotations\', files(ii).name]);
    
    mask = roipoly(mass.mass_ROI, mass.mass_outline(:,1), mass.mass_outline(:,2));
    for jj = 1:49; mask = imdilate(mask, strel('disk', 1)); end
    
    p_list = find(~mask);
    m_sm = imfilter(mass.mass_ROI, fspecial('gaussian', 90, 18), 'symmetric');
    m_su = m_sm - mass.mass_spline_it;
    m_su(p_list) = 0;
    outline_x = mass.mass_outline(:,1);
    outline_y = mass.mass_outline(:,2);
    
    figure('WindowStyle', 'docked', 'Name', mass.name);
    subplot(1,2,1)
    image(mass.mass_ROI - m_su); colormap(gray(256)); axis('image');
    hold('on');
    %plot(outline_x([1:end, 1]), outline_y([1:end, 1]), 'k:');
    subplot(1,2,2)
    image(mass.mass_ROI - uint8(mass.mass_sub_smooth)); colormap(gray(256)); axis('image');
    %imagesc(double(mass.mass_spline) - double(mass.mass_spline_it)); axis('image');
    hold('on');
    %plot(outline_x([1:end, 1]), outline_y([1:end, 1]), 'k:');
    colormap(jet(256));
    clear m_sm m_su mass;
end
%%
for ii = 161:185
    load(['..\..\dev\annotations\', files(ii).name]);
    mask = roipoly(mass.mass_ROI, mass.mass_outline(:,1), mass.mass_outline(:,2));
    for jj = 1:49; mask = imdilate(mask, strel('disk', 1)); end
    
    p_list = find(~mask);
    m_sm = imfilter(mass.mass_ROI, fspecial('gaussian', 90, 18), 'symmetric');
    m_su = m_sm - mass.mass_spline_it;
    m_su(p_list) = 0;
    try
        mbg = bg.dense(1:size(m_su, 1), 1:size(m_su, 2)) + m_su;
        figure('WindowStyle', 'docked', 'Name', mass.name);
    end
    image(mbg); colormap(gray(256)); axis('image');
end
%%
for ii = 12:185
    load(['..\..\dev\annotations\', files(ii).name]);
    mass.mass_spline_it = [];
    mass = rmfield(mass, {'mass_list20', 'mass_list20_40'});
    save(['..\..\dev\annotations\', files(ii).name], 'mass');
end