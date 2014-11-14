% run these commands after x_b, x_sw have been calculated to recalculate x_g with different x_b.
% also creates an image superimposing the x_g<0 region onto the original image

%b_fac=1.0;
disp('x_b scaled by (%):');
disp(b_fac*100);
x_g_2 = (x_b > 0).*(x_sw - CAL_A(wedgesizes(i_file),ikv)*x_b*b_fac - CAL_B(wedgesizes(i_file),ikv))./(CAL_D(wedgesizes(i_file),ikv)*x_b*b_fac + CAL_E(wedgesizes(i_file),ikv)) ;
imfinal = uint8(((double(IMAGE)/255 - double(x_g_2<0))>0).*(double(IMAGE)/255 - double(x_g_2<0))*255);
imshow(imfinal,map)    

% work out area and volume of breast
area_b = sum(sum(x_b > 0)); % nb. this is in units of (pixels)^2
area_b_nm = sum(sum(x_b.*(1-marker_mask) > 0));
area_g = sum(sum(x_g_2 > 0));
area_g_nm = sum(sum(x_g_2.*(1-marker_mask) > 0));

area_gt_5 = sum(sum(x_g_2.*(1-marker_mask) > 5 ));

volume_b = sum(sum(x_b))*pixel_area*1e-3*area_resize_factor; % nb. units are m^3 
volume_b_nm = sum(sum(x_b.*(1-marker_mask)))*pixel_area*1e-3*area_resize_factor;
volume_g = sum(sum( (x_g_2 > 0).*x_g ))*pixel_area*1e-3*area_resize_factor;
volume_g_nm = sum(sum( (x_g_2.*(1-marker_mask) > 0).*x_g_2 ))*pixel_area*1e-3*area_resize_factor;
    
dense_by_area_nm = 100 * area_g_nm / area_b_nm ;
dense_by_vol_nm  = 100 * volume_g_nm / volume_b_nm ;

dense_by_area_gt_5 = 100 * area_gt_5 / area_b_nm;

dense_by_area = 100 * area_g / area_b ;
dense_by_vol = 100 * volume_g / volume_b ;
fat_by_vol_nm=100.-dense_by_vol_nm;
volume_f_nm=volume_b-volume_g_nm;    

%disp('dense area(%)   area_gt_5 (%)  dense vol(%)   vol_g (cm^3)   vol_b (cm^3)  ');
disp('Vol_b    vol_g    vol_f    dense vol(%) fat vol(%) dense area(%)   area_gt_5 (%)');
disp([volume_b*100*100*100 volume_g_nm*100*100*100 volume_f_nm*100*100*100 dense_by_vol_nm fat_by_vol_nm dense_by_area_nm dense_by_area_gt_5 ]);
    
