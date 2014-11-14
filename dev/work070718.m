for ii = 1:183
    
    load(['C:\isbe\dev\annotations\', files(ii).name]);
    f_name = ['m', files(ii).name(6:end-3), 'jpg'];
    imwrite(mass.mass_ROI,...
        ['C:\isbe\dev\images\ROI\', f_name],...
        'quality', 100);
    imwrite(mass.mass_sub_smooth,...
        ['C:\isbe\dev\images\mass\', f_name],...
        'quality', 100);
    
    bg = uint8(double(mass.mass_ROI) - mass.mass_sub_smooth);
    imwrite(bg,...
        ['C:\isbe\dev\images\bg\', f_name],...
        'quality', 100);
    clear mass bg;
end
%%
for ii = 1:183
    
    f_name = ['m', files(ii).name(6:end-3), 'jpg'];
    b_new = ['m', files(ii).name(6:end-4), '_b.jpg'];
    m_new = ['m', files(ii).name(6:end-4), '_m.jpg'];
    r_new = ['m', files(ii).name(6:end-4), '_r.jpg'];
    movefile(['C:\isbe\dev\images\bg\', f_name],...
        ['C:\isbe\dev\images\', b_new]);
    movefile(['C:\isbe\dev\images\mass\', f_name],...
        ['C:\isbe\dev\images\', m_new]);
    movefile(['C:\isbe\dev\images\ROI\\', f_name],...
        ['C:\isbe\dev\images\', r_new]);
end
%%
for ii = 1:183
    
    m_new = ['m', files(ii).name(6:end-4), '_m.jpg'];
    b_new = ['m', files(ii).name(6:end-4), '_b.jpg'];
    movefile(['C:\isbe\dev\images\bg\', m_new],...
        ['C:\isbe\dev\images\', b_new]);
end
%%
for ii = 1:179
    
    load(['C:\isbe\dev\annotations\', files(ii).name]);
    f_name = ['m', files(ii).name(6:end-4), '_b.jpg'];
    
    mass.mass_sub_it = mass.mass_sub_smooth + double(mass.mass_spline)...
        - double(mass.mass_spline_it);
    save(['C:\isbe\dev\annotations\', files(ii).name],...
        'mass');
    
    bg = uint8(double(mass.mass_ROI) - mass.mass_sub_it);
    
    imwrite(bg,...
        ['C:\isbe\dev\images\', f_name],...
        'quality', 100);
    clear mass bg;
end
%
%%
for ii = 1:183
    display([files(ii).name(6:end-4)]);
    display([' INCLUDEPICTURE  "images/m', files(ii).name(6:end-4), '_r.jpg" \d  \* MERGEFORMAT ']);
    display([' INCLUDEPICTURE  "images/m', files(ii).name(6:end-4), '_b.jpg" \d  \* MERGEFORMAT ']);
end
%%
for ii = 1:183
    load(['C:\isbe\dev\annotations\', files(ii).name]);
    
    i1 = imread(mass.name);
    [rows cols] = size(i1); clear i1;
    if (mass.R1 == 1) || (mass.C1 == 1) || (mass.R2 == rows) || (mass.C2 == cols)
        display([num2str(ii)]);
        %figure('WindowStyle', 'docked', 'name', files(ii).name);
        %imagesc(mass.mass_ROI); colormap(gray(256)); axis image;
    end
    clear mass;
end
%%
for ii = 1:183
    load(['C:\isbe\dev\annotations\', files(ii).name]);
    mass.name = ['C:\isbe\mammograms\new_CAD\BMP_2004\o', files(ii).name(3:11), '.bmp'];
    save(['C:\isbe\dev\annotations\', files(ii).name], 'mass');
end
%%