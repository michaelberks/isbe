for ii = 0:5:25
    
    imwrite((p_thickness > ii/100), ['gland_area', zerostr(ii,2), '.jpg']);
end
write_im_from_colormap(gland_thickness, 'gland_thickness.jpg', gray(256), [0 30]); 