% Stepwedge calibration data
% Use this to generate a relationship between gland thickness, breast
% thickness and stepwedge thickness

% x-coordinates = breast thickness

xb1 = [20 20 20 20 20];
xb2 = [30 30 30 30 30 30];
xb3 = [40 40 40 40 40 40];
xb4 = [50 50 50 50 50 50 50 50 50];
xb5 = [60 60 60 60 60 60 60 60];
xb6 = [70 70 70 70 70 70 70 70 70 70];

% y-coordinates = gland thickness

xg1 = [20 15 10 5 0];
xg2 = [30 25 20 15 10 0];
xg3 = [40 35 30 25 20 10];
xg4 = [50 40 35 30 25 20 15 10 0];
xg5 = [50 45 40 35 20 15 10 0];
xg6 = [45 40 35 30 25 20 15 10 5 0];

% z-coordinates = stepwedge thickness

sw1 = [1.51 1.34 1.14 0.97 0.81];
sw2 = [2.50 2.30 2.06 1.87 1.67 1.29];
sw3 = [3.81 3.67 3.31 3.01 2.75 2.24];
sw4 = [6.86 5.21 4.82 4.21 3.90 3.71 3.33 2.95 2.43];
sw5 = [9.07 7.59 6.70 6.16 4.60 4.18 3.78 3.29];
sw6 = [13.99 10.60 8.63 7.54 6.66 6.02 5.35 4.83 4.38 3.97];

figure(1);
plot(xb1,sw1,'bo')
title('Breast thickness versus stepwedge thickness')
xlabel('Breast thickness')
ylabel('Stepwedge thickness')
hold on;
plot(xb2,sw2,'ro')
plot(xb3,sw3,'co')
plot(xb4,sw4,'go')
plot(xb5,sw5,'mo')
plot(xb6,sw6,'ko')

figure(2);
plot(xg1, sw1,'bo')
title('Glandular thickness versus stepwedge thickness')
xlabel('Glandular thickness')
ylabel('Stepwedge thickness')
axis( [-5, Inf, -Inf, Inf] )
hold on;
plot(xg2, sw2,'ro')
plot(xg3, sw3,'co')
plot(xg4, sw4,'go')
plot(xg5, sw5,'mo')
plot(xg6, sw6,'ko')

figure(3);
plot3(xb1, xg1, sw1,'bo')
title('Relationship between stepwedge thickness, glandular thickness and total breast thickness')
xlabel('Breast thickness')
ylabel('Glandular thickness')
zlabel('Stepwedge thickness')
hold on;
plot3(xb2, xg2, sw2,'ro')
plot3(xb3, xg3, sw3,'co')
plot3(xb4, xg4, sw4,'go')
plot3(xb5, xg5, sw5,'mo')
plot3(xb6, xg6, sw6,'ko')

xb = 20:10:70;
xg = 0:5:50;
sw = [0.81 1.29 nan 2.43 3.29 3.97; 0.97 nan nan nan nan 4.38; ...
    1.14 1.67 2.24 2.95 3.78 4.83; 1.34 1.87 nan 3.33 4.18 5.35; ...
    1.51 2.06 2.75 3.71 4.60 6.02; nan 2.30 3.01 3.90 nan 6.66; ...
    nan 2.50 3.31 4.21 nan 7.54; nan nan 3.67 4.82 6.16 8.63; ...
    nan nan 3.81 5.21 6.70 10.60; nan nan nan nan 7.59 13.99; ...
    nan nan nan 6.86 9.07 nan];

figure(4);
[x y] = meshgrid(xb, xg);
mesh(xb, xg, sw)
title('Relationship between stepwedge thickness, glandular thickness and total breast thickness')
xlabel('Breast thickness')
ylabel('Glandular thickness')
zlabel('Stepwedge thickness')
axis( [-Inf, Inf, -5, Inf, -Inf, Inf] )

figure(5);
[xi yi] = meshgrid(20:1:70, 0:1:50);
zi = interp2(xb, xg, sw, xi, yi);
mesh(xi, yi, zi)
axis( [-Inf, Inf, -1, Inf, -Inf, Inf] )

% is this the right thing to do?
% Ideally we require xg = f(xb, sw) rather than sw = f(xb, xg)
% try non-linear interpolation
% even if you can interpolate, you still need to generate the equation for
% the function