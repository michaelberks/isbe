% Stepwedge calibration data
% Use this to generate a relationship between gland thickness, breast
% thickness and stepwedge thickness

% Start by interpolating sw vs xg

% % x-coordinates = gland thickness
% 
% xg20 = [0 5 10 15 20];
% xg30 = [0 10 15 20 25 30];
% xg40 = [10 20 25 30 35];
% xg50 = [0 10 15 20 25 30 35 40 50];
% xg60 = [0 10 15 20 35 40 45 50];
% xg70 = [0 5 10 15 20 25 30 35 40 45];
% 
% % y-coordinates = stepwedge thickness
% 
% sw20 = [0.81 0.97 1.14 1.34 1.51];
% sw30 = [1.29 1.67 1.87 2.06 2.30 2.50];
% sw40 = [2.24 2.75 3.01 3.31 3.67];
% sw50 = [2.43 2.95 3.33 3.71 3.90 4.21 4.82 5.21 6.86];
% sw60 = [3.29 3.78 4.18 4.60 6.16 6.70 7.59 9.07];
% sw70 = [3.97 4.38 4.83 5.35 6.02 6.66 7.54 8.63 10.60 13.99];

% figure(1);
% plot(xg20, sw20,'bo')
% title('Glandular thickness versus stepwedge thickness')
% xlabel('Glandular thickness')
% ylabel('Stepwedge thickness')
% axis( [-5, Inf, -Inf, Inf] )
% hold on;
% plot(xg30, sw30,'ro')
% plot(xg40, sw40,'co')
% plot(xg50, sw50,'go')
% plot(xg60, sw60,'mo')
% plot(xg70, sw70,'ko')
% 
% xg_int = 0:50;
% sw20_int = interp1(xg20, sw20, xg_int,'linear','extrap');
% sw30_int = interp1(xg30, sw30, xg_int,'linear','extrap');
% sw40_int = interp1(xg40, sw40, xg_int,'spline','extrap');
% sw50_int = interp1(xg50, sw50, xg_int,'spline','extrap');
% sw60_int = interp1(xg60, sw60, xg_int,'spline','extrap');
% sw70_int = interp1(xg70, sw70, xg_int,'spline','extrap');
% 
% % view the interpolated fits
% % all have been extrapolated to have points between 0 and 50mm gland
% plot(xg_int, sw20_int)
% plot(xg_int, sw30_int)
% plot(xg_int, sw40_int)
% plot(xg_int, sw50_int)
% plot(xg_int, sw60_int)
% plot(xg_int, sw70_int)
% 
% % now use the interpolated sw values and view as a mesh
% xb = 20:10:70;
% sw_int = [sw20_int; sw30_int; sw40_int; sw50_int; sw60_int; sw70_int];
% 
% figure(2);
% [x y] = meshgrid(xg_int,xb);
% mesh(xg_int, xb, sw_int)
% title('Relationship between stepwedge thickness, glandular thickness and total breast thickness')
% xlabel('Glandular thickness')
% ylabel('Breast thickness')
% zlabel('Stepwedge thickness')
% 
% 
% figure(3);
% [xi yi] = meshgrid(0:0.5:50, 20:0.5:70);
% zi = interp2(xg_int, xb, sw_int, xi, yi);
% mesh(xi, yi, zi)
% title('Relationship between stepwedge thickness, glandular thickness and total breast thickness')
% xlabel('Glandular thickness (mm)')
% ylabel('Breast thickness (mm)')
% zlabel('Stepwedge thickness (mm)')

% This mesh has been generated for interest but cannot be used in the main
% programme.
% The main 'stepwedge' programme requires a glandular thickness to be 
% generated at every pixel based on the total thickness and stepwedge
% thickness.
% These are both known following the marker and stepwedge mark-up
%%


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
gt = zeros(6, length(step_height_pts));
gt(1,:) = interp1(sw20, xg20, step_height_pts,'linear','extrap'); 
gt(2,:) = interp1(sw30, xg30, step_height_pts,'linear','extrap');
gt(3,:) = interp1(sw40, xg40, step_height_pts,'linear','extrap');
gt(4,:) = interp1(sw50, xg50, step_height_pts,'linear','extrap');
gt(5,:) = interp1(sw60, xg60, step_height_pts,'linear','extrap');
gt(6,:) = interp1(sw70, xg70, step_height_pts,'linear','extrap');


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
    0.2:0.2:14, (20:70)', 'linear');


figure;
mesh(0.2:0.2:14, (20:70)', gland_thickness);
axis([0 14 20 70 0 70]);
%%
figure; hold on;
plot(step_height_pts, gt(20,:), 'r-x');
plot(step_height_pts, gt(30,:), 'g-x');
plot(step_height_pts, gt(40,:), 'b-x');
plot(step_height_pts, gt(50,:), 'm-x');
plot(step_height_pts, gt(60,:), 'c-x');
plot(step_height_pts, gt(70,:), 'y-x');

figure; hold on;
plot(sw20, 20, 'x');
plot(sw30, 30, 'x');
plot(sw40, 40, 'x');
plot(sw50, 50, 'x');
plot(sw60, 60, 'x');
plot(sw70, 70, 'x');

