function [x_g,x_sw] = sw_transform(IMAGE,x_b,sw_lookup, CAL_A,CAL_B,CAL_D,CAL_E)

figure(4);
plot(sw_lookup);

% now transform IMAGE to x_sw within breast area, and zero outside

% need a version of the image with pixel values from 1 to 256,
% so these values can be used to address the sw_lookup array.
tempIMAGE = double(IMAGE)+1;

x_sw = sw_lookup(tempIMAGE).* (x_b > 0);

clear tempIMAGE;

figure(6);
display_array(x_sw);

% now transform original image (within breast area) with 
% x_g = (x_sw - CAL_A*x_b - CAL_B) / (CAL_D*x_b + CAL_E)

x_g = (x_b > 0).*(x_sw - CAL_A*x_b - CAL_B)./(CAL_D*x_b + CAL_E) ;
