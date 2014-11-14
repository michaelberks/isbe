% Stepwedge calibration data
% Use this to generate a relationship between gland thickness, breast
% thickness and stepwedge thickness

function [gland_thickness] = sw_vs_gland(sw_lookup, breast_thickness)

%To compute for multiple values make sure input arguments are NOT
%row/column otherwise will calculate mesh
% x-coordinates = gland thickness

xg20 = [0 5 10 15 20];
xg30 = [0 10 15 20 25 30];
xg40 = [10 20 25 30 35];
xg50 = [0 10 15 20 25 30 35 40 50];
xg60 = [0 10 15 20 35 40 45 50];
xg70 = [0 5 10 15 20 25 30 35 40 45];

% y-coordinates = stepwedge thickness

sw20 = [0.81 0.97 1.14 1.34 1.51];
sw30 = [1.29 1.67 1.87 2.06 2.30 2.50];
sw40 = [2.24 2.75 3.01 3.31 3.67];
sw50 = [2.43 2.95 3.33 3.71 3.90 4.21 4.82 5.21 6.86];
sw60 = [3.29 3.78 4.18 4.60 6.16 6.70 7.59 9.07];
sw70 = [3.97 4.38 4.83 5.35 6.02 6.66 7.54 8.63 10.60 13.99];

step_height_pts = 0.5:0.5:14;
gt(1,:) = interp1(sw20, xg20, step_height_pts,'linear','extrap'); 
gt(2,:) = interp1(sw30, xg30, step_height_pts,'linear','extrap');
gt(3,:) = interp1(sw40, xg40, step_height_pts,'spline','extrap');
gt(4,:) = interp1(sw50, xg50, step_height_pts,'spline','extrap');
gt(5,:) = interp1(sw60, xg60, step_height_pts,'spline','extrap');
gt(6,:) = interp1(sw70, xg70, step_height_pts,'spline','extrap');


% do this for all stepheights e.g. (gt1 > 20) = 20, (gt2 > 30) = 30;
gt(1, gt(1,:) > 20) = 20;
gt(2, gt(2,:) > 30) = 30;
gt(3, gt(3,:) > 40) = 40;
gt(4, gt(4,:) > 50) = 50;
gt(5, gt(5,:) > 60) = 60;
gt(6, gt(6,:) > 70) = 70;
gt(gt < 0) = 0;


gland_thickness = interp2(repmat(step_height_pts, 6, 1),...
    repmat((20:10:70)', 1, 28),...
    [gt(1,:); gt(2,:); gt(3,:); gt(4,:); gt(5,:); gt(6,:)],...
    sw_lookup, breast_thickness, 'linear');
