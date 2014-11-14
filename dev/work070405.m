%%
for ii=1:16
    mass(ii) = eval(['mass', zerostr(ii,3)]);
end
%%
mass016.name = 'C:\isbe\mammograms\new_CAD\BMP_2004\o04_125LCC.bmp';
mass016.target_name = 'C:\isbe\mammograms\new_CAD\BMP_2004\o04_125RCC.bmp';
mass016.target_ROI = target_ROI;
mass016.target_pos = [1500, 2150, 3280, 3944];
mass016
%%
% %%
% load C:\isbe\dev\mass_data
% mass1 = mass_data(1);
% clear mass_data
% mammo1 = fliplr(imread('C:\isbe\mammograms\new_CAD\BMP_2004\o04_004LML.bmp'));
% mammo1 = double(mammo1) ./ 255;
% no_mass1 = mass1.mass_ROI;
% mask1 = poly2mask(mass1.mass_outline(:,1), mass1.mass_outline(:,2), size(mass1.mass_ROI, 1), size(mass1.mass_ROI,2));
% no_mass1(mask1) = NaN;
% no_mass1(654:807,:) = NaN;
% no_mass1 = double(no_mass1) ./ 255;
% weight1 = fspecial('gaussian', 807, 20);
% figure; imagesc(mammo1); colormap(gray); axis image;
% %%
% mammo2 = mammo1(500:end, 500:end); clear mammo1;
% mammo1 = mammo2; clear mammo2;