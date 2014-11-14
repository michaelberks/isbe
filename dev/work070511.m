%%
for ii = 1:15
    load(sub_masses(ii).name);
    figure('WindowStyle', 'docked');
    temp = bg.dense;
    temp(mass.p_list1500) = temp(mass.p_list1500) + uint8(mass.mass_sub);
    image(temp); colormap(gray(256)); axis image; hold on; clear temp;
    figure('WindowStyle', 'docked');
    temp = bg.fatty;
    temp(mass.p_list1500) = temp(mass.p_list1500) + uint8(mass.mass_sub);
    image(temp); colormap(gray(256)); axis image; hold on; clear temp;
    figure('WindowStyle', 'docked');
    temp = bg.parenchyma;
    temp(mass.p_list1500) = temp(mass.p_list1500) + uint8(mass.mass_sub);
    image(temp); colormap(gray(256)); axis image; hold on; clear temp;
    figure('WindowStyle', 'docked');
    temp = bg.streaky;
    temp(mass.p_list1500) = temp(mass.p_list1500) + uint8(mass.mass_sub);
    image(temp); colormap(gray(256)); axis image; hold on; clear temp;
    %plot(mass.mass_outline([1:end,1],1)-mass.offset1500(1)+750,...
    %    mass.mass_outline([1:end,1],2)-mass.offset1500(2)+750, 'r');
    clear mass;
end
%%
    %figure('WindowStyle', 'docked');
    %load(sub_masses(2).name);
    subplot(2,4,1)
    temp = bg.dense;
    image(temp); colormap(gray(256)); axis image; hold on; clear temp;
    subplot(2,4,2)
    temp = bg.fatty;
    image(temp); colormap(gray(256)); axis image; hold on; clear temp;
    subplot(2,4,3)
    temp = bg.parenchyma;
    image(temp); colormap(gray(256)); axis image; hold on; clear temp;
    subplot(2,4,4)
    temp = bg.streaky;
    image(temp); colormap(gray(256)); axis image; hold on; clear temp;
%%
cd C:\isbe\dev\subtraction
for ii = 1:15
    load(sub_masses(ii).name);
    mass.offset1500 = round(mass.offset1500);
    mass_bw = poly2mask(mass.mass_outline(:,1)-mass.offset1500(1)+750,...
                        mass.mass_outline(:,2)-mass.offset1500(2)+750,...
                        1500, 1500);
    for jj = 1:20
        mass_bw = imdilate(mass_bw, strel('disk', 1));
    end
    
    mass_rp = regionprops(bwlabel(mass_bw,4), 'PixelIdxList');
    mass.p_list1500 = mass_rp.PixelIdxList;
    clear mass_rp mass_bw;
    save(sub_masses(ii).name, 'mass');
    figure('WindowStyle', 'docked');
    try
        temp = zeros(1500,1500);
        temp(mass.p_list1500) = uint8(mass.mass_sub);
        image(temp); colormap(gray(256)); axis image;
    catch
        display(sub_masses(ii).name);
    end
 
    clear temp mass;
end
%%
cd C:\isbe\dev\subtraction\annotations\
figure('WindowStyle', 'docked');
for ii = 1:15
    load(sub_masses(ii).name);
    mass_bw = poly2mask(mass.mass_outline(:,1)-mass.offset1500(1)+750,...
                        mass.mass_outline(:,2)-mass.offset1500(2)+750,...
                        1500, 1500);
    subplot(3,5,ii);
    imagesc(mass_bw); colormap gray; axis image;
 
    clear mass mass_bw;
end