L = 246;
base_l = 280;
d = zeros(100,1);
d(1) = base_l;

theta_min = pi/12;

theta = linspace(theta_min, pi/2, 50);
theta = [fliplr(theta) theta];

xc = d(1);% - k / tan(theta(51));
yc = 0;

xa = 0;
ya = L;
xb = zeros(100,1);
yb = zeros(100,1);

f1 = figure('Position', [100,100, 600, 400], 'WindowStyle', 'normal', 'Color', [1, 1, 1]);
hold on; 
axis equal off; axis([-5 600 -5 400]);
%title('Camera arm system');

plot([0 0], [0 ya], 'k', 'LineWidth', 2.0);
plot([0 2*base_l], [0 0], 'k');
%text(-50, 50, {'Patient'; 'here'});

for i_t = 1:100
    
    %xa(i_t) = 0; by definition
    
    xb(i_t) = xc + L*cos(theta(i_t)); 
    yb(i_t) = yc + L*sin(theta(i_t));
    
    d(i_t) = sqrt(xb(i_t)^2 + (yb(i_t)-ya).^2);
    
    base_x = [xa xa + base_l*(xb(i_t) - xa)/d(i_t)];
    base_y = [ya ya + base_l*(yb(i_t) - ya)/d(i_t)];
    acuator_x = [xa xb(i_t)];
    acuator_y = [ya yb(i_t)];
    
    camera_x = [xc xc + 2*(xb(i_t) - xc)];
    camera_y = [yc yc + 2*(yb(i_t) - yc)];
    
    if i_t == 1

        plot(acuator_x, acuator_y, 'm', 'LineWidth', 2.0,...
            'XDataSource', 'acuator_x',...
            'YDataSource', 'acuator_y');
        plot(base_x, base_y, 'k', 'LineWidth', 2.0,...
            'XDataSource', 'base_x',...
            'YDataSource', 'base_y');
        
        plot(camera_x, camera_y, 'b', 'LineWidth', 2.0,...
            'XDataSource', 'camera_x',...
            'YDataSource', 'camera_y');
        
        plot(xb(i_t), yb(i_t), 'g.', 'MarkerSize', 8,...
            'XDataSource', 'xb(i_t)',...
            'YDataSource', 'yb(i_t)');
        
        plot(xa, ya, 'r.', 'MarkerSize', 8);
    else
        refreshdata(gca, 'caller');
    end
    
    frame1 = getframe(f1);
    gif1 = frame2im(frame1);
    [gif1a map] = rgb2ind(gif1, 2^16);
    
    if i_t == 1
        imwrite(gif1a, map, 'camera_arm.gif', 'gif', 'WriteMode', 'overwrite', 'Loop', Inf, 'DelayTime', 0.05);
    else
        imwrite(gif1a, map, 'camera_arm.gif', 'gif', 'WriteMode', 'append', 'DelayTime', 0.05);
    end
    
end
%%
%Computing the forces involved
F = 200; %Max actuator force in N
M = 2; %Mass of camera in kg
g = 10; %Gravity
L1 = L; %Length of the post to the attachment of the actuator
L2 = 300; %Length of the post to the point mass of the cameras

phi = theta_min + atan((ya - yb(50)) / xb(50));

Mc = M*g*cos(theta_min)*L2;
Ma = F*sin(phi)*L1;

max_weight = Ma / (g*cos(theta_min)*L2);


