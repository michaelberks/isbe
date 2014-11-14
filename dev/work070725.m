mean_mass = zeros(179,1);
for ii = 1:179
    load(['C:\isbe\dev\annotations\', files(ii).name]);
    mean_mass(ii) = sum(mass.mass_sub_it(:)) / length(mass.mass_list);
    clear mass;
end
%%
for ii = 161:179
    load(['C:\isbe\dev\annotations\', files(ii).name]);
    figure('WindowStyle', 'docked')
    imagesc(uint8(mass.mass_sub)); axis('image'); colormap(gray(256));
    hold on;
    plot(mass.mass_outline(:,1), mass.mass_outline(:,2), 'r:');
end
%%
gf1 = fspecial('gaussian', 10, 0.5);
gf2 = fspecial('gaussian', 10, 1);
for ii = 10:20
    load(['C:\isbe\dev\annotations\', files(ii).name]);
    
    ms1 = imfilter(mass.mass_ROI, gf1, 'symmetric');
    ms2 = imfilter(mass.mass_ROI, gf2, 'symmetric');
    dm = double(ms2) - double(ms1);
    
    mass_sub = mass.mass_sub_it + dm;
    mass_sub(setdiff(1:end, mass.mass_list)) = 0;
    
    figure('WindowStyle', 'docked');
    subplot(1,2,1)
    imagesc(dm); axis('image'); colormap(gray(256));
    subplot(1,2,2)
    imagesc(uint8(mass_sub)); axis('image'); colormap(gray(256));
    clear mass;
end
%%
for ii = 1:20
    load(['C:\isbe\dev\annotations\', files(ii).name]);
    
    [r c] = size(mass.mass_ROI)
    if r < 1500 && c < 1500
        nm = bg.parenchyma(1:r, 1:c) + uint8(mass.mass_sub_it);
        figure('WindowStyle', 'docked'); imagesc(nm); axis image; colormap(gray(256));
    end
end


    

