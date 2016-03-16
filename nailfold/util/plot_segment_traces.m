function plot_segment_traces(pixels_per_mm, segments_s, segments_ns, motor_x, motor_y, motor_z, sharpness)

xy_speed = [0; sqrt(diff(motor_x).^2 + diff(motor_y).^2)];
z_speed = [0; diff(motor_z)];   
frame_nums = 1:length(xy_speed);

figure;
a1 = subplot(3,1,1);hold on; 
plot(xy_speed, 'b');
plot(abs(z_speed), 'r');
xlabel('Time', 'fontsize', 18);
ylabel('Motor speed', 'fontsize', 18);
legend({'|XY| speed', '|Z| speed'});
set(a1, 'ylim', [ 0 0.04]);

%plot(sharpness_change / (max(sharpness_change)-min(sharpness_change)), 'g--');
a2 = subplot(3,1,2); hold on;
xlabel('Time', 'fontsize', 18);
ylabel('Z position', 'fontsize', 18);
plot(0,0,'r-');
plot(0,0,'b-');
legend(a2, {'Stationary/moving in Z', '|XY| movement'});

% a3 = subplot(3,1,3); hold on;
% xlabel('Time');
% ylabel('Sharpness');
% plot(0,0,'r-');
% plot(0,0,'b-');
% legend(a3, {'Stationary/moving in Z', 'XY movement'});

for i_sg = 1:length(segments_s)
    plot(a2, frame_nums(segments_s{i_sg}), motor_z(segments_s{i_sg}), 'r.');
    %plot(a3, frame_nums(segments_s{i_sg}), sharpness(segments_s{i_sg}), 'r.');
end

for i_sg = 1:length(segments_ns)
    plot(a2, frame_nums(segments_ns{i_sg}), motor_z(segments_ns{i_sg}), 'b.');
    %plot(a3, frame_nums(segments_ns{i_sg}), sharpness(segments_ns{i_sg}), 'b.');
end
set(a2, 'ylim', [ min(motor_z) max(motor_z)]);

hw = 320 / pixels_per_mm;
hh = 240 / pixels_per_mm;
frame_rect_x = [-hw hw hw -hw -hw];
frame_rect_y = [-hh -hh hh hh -hh];

%figure; hold all; axis equal; a1 = gca;

a3 = subplot(3,1,3); hold all; axis equal;
plot(0,0,'b-');
plot(0,0,'r-');
plot(0,0,'rx');
legend(a3, {'Motor path', 'Stationary frames', 'Stationary centres'});
plot(a3, motor_x, motor_y, 'b.', 'markersize', 4);
%
min_x = inf;
max_x = -inf;
min_y = inf;
max_y = -inf;
for i_seg = 1:length(segments_s)

    frame_idx_i = segments_s{i_seg};

    seg_x = motor_x(frame_idx_i);
    seg_y = motor_y(frame_idx_i);

    cx = median(seg_x);
    cy = median(seg_y);

    plot(a3, cx, cy, 'rx');
    plot(a3, cx+frame_rect_x, cy+frame_rect_y, 'r', 'linewidth', 2); 
    min_x = min(min_x,min(cx+frame_rect_x));
    max_x = max(max_x,max(cx+frame_rect_x));
    min_y = min(min_y,min(cy+frame_rect_y));
    max_y = max(max_y,max(cy+frame_rect_y));   
    
end
axis([min_x max_x min_y max_y]);

xlabel('X position', 'fontsize', 18);
ylabel('Y position', 'fontsize', 18);
set(a1, 'fontsize', 18);
set(a2, 'fontsize', 18);
set(a3, 'fontsize', 18);
%%
figure('windowstyle', 'normal');
hold all; axis equal; a3 = gca;

plot(0,0,'b-.');
plot(0,0,'r-');
plot(0,0,'rx');
legend(a3, {'Motor path', 'Stationary segment frames', 'Stationary segment centres'}, 'location', 'north');
plot(a3, motor_x, motor_y, 'b.', 'markersize', 8);
%
min_x = inf;
max_x = -inf;
min_y = inf;
max_y = -inf;
for i_seg = 1:length(segments_s)

    frame_idx_i = segments_s{i_seg};

    seg_x = motor_x(frame_idx_i);
    seg_y = motor_y(frame_idx_i);

    cx = median(seg_x);
    cy = median(seg_y);

    plot(a3, cx, cy, 'rx');
    plot(a3, cx+frame_rect_x, cy+frame_rect_y, 'r', 'linewidth', 2); 
    min_x = min(min_x,min(cx+frame_rect_x));
    max_x = max(max_x,max(cx+frame_rect_x));
    min_y = min(min_y,min(cy+frame_rect_y));
    max_y = max(max_y,max(cy+frame_rect_y));   
    
end
axis([min_x max_x min_y max_y]);

xlabel('X position', 'fontsize', 14);
ylabel('Y position', 'fontsize', 14);
set(a3, 'fontsize', 14);