function h=plot_circle(green, a, b, c)

pts = linspace(0, 2*pi, 200);
source_x = a*cos(pts) + c; %source landmark points
source_y = a*sin(pts) + c;

int_x = b*cos(pts) + c; %points to interpolate
int_y = b*sin(pts) + c;

temp1 = [repmat(a, 1, 50), linspace(a, -a, 50), repmat(-a, 1, 50), linspace(-a, a, 50)];
temp2 = [linspace(-a, a, 50), repmat(a, 1, 50), linspace(a, -a, 50), repmat(-a, 1, 50)];

target_x = [temp1(26:200), temp1(1:25)] + c; %target landmark points
target_y = [temp2(26:200), temp2(1:25)] + c;

%translate source x and y about 0
source_x = source_x - mean(source_x);
source_y = source_y - mean(source_y);

%translate target x and y about 0
target_x = target_x - mean(target_x);
target_y = target_y - mean(target_y);

%translate int x and y about 0
mx = mean(int_x); my = mean(int_y);
int_x = int_x - mx;
int_y = int_y - my;

%scale source, int, target x and y to fit inside unit circle
scale_src = sqrt(max(source_x.^2 + source_y.^2));
scale_target = sqrt(max(target_x.^2 + target_y.^2));
scale_factor = max(scale_src, scale_target)*2;

source_x = source_x / scale_factor;
source_y = source_y / scale_factor;

target_x = target_x / scale_factor;
target_y = target_y / scale_factor;

int_x = int_x / scale_factor;
int_y = int_y / scale_factor;

L_inv = tps_weights(source_x, source_y, green);
new_x = tps_warp(L_inv, source_x, source_y, target_x, int_x, int_y, green);
new_y = tps_warp(L_inv, source_x, source_y, target_y, int_x, int_y, green);

h = [new_x', new_y'];

figure;
hold on;

plot(target_x, target_y, 'b-');  
plot(source_x, source_y, 'r-');
plot(int_x, int_y, 'bx');
plot(new_x, new_y, 'rx');