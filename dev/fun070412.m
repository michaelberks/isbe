function fun070412(data_in, outpath)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make cut and paste mass images
%

N = length(data_in);
gauss_filt = fspecial('gaussian', 25, 5);

for ii = 1:N
    mass = data_in(ii);
    
    new_mass_ROI1 = mass.target_ROI + uint8(fliplr(mass.mass_sub));
    new_mass_ROI2 = mass.target_ROI + ...
        imfilter(uint8(fliplr(mass.mass_sub)), gauss_filt);
    new_mass_ROI3 = uint8(double(mass.target_spline)...
        + fliplr(mass.mass_sub)); 
    
    mammo = imread(mass.target_name);
    r1 = mass.target_pos(1); r2 = mass.target_pos(2);
    c1 = mass.target_pos(3); c2 = mass.target_pos(4);
    
    mammo(r1:r2, c1:c2) = new_mass_ROI1;
    imwrite(mammo, [outpath, 'unsmoothed', zerostr(ii, 3), '.bmp']);
    
    mammo(r1:r2, c1:c2) = new_mass_ROI2;
    imwrite(mammo, [outpath, 'smoothed', zerostr(ii, 3), '.bmp']);
    
    mammo(r1:r2, c1:c2) = new_mass_ROI3;
    imwrite(mammo, [outpath, 'splined', zerostr(ii, 3), '.bmp']);
end