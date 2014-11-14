%Fig 1: Build 1D step-edges
step_edges1D = 2*ones(64, 64);
for offset = 1:64;
    step_edges1D(1:offset, offset) = 1;
end
%
%Construct 1D dual-tree for each offset step-edge
[dummy dt] = dtwavexfm(step_edges1D, 4, 'near_sym_b','qshift_b');
%
%Look at phase of level 1 coefficients - should be linear in support range
% and zero at edge
figure;
for coeff = 1:8
    subplot(4, 2, coeff); hold on; xlabel(['Coeff: ', num2str(coeff)]);
    plot(1:64, angle(dt{3}(coeff,:)), 'b');
    plot(1:64, angle(dt{3}(coeff,:)), 'rx');
    plot([1 64], [0 0], 'r-');
end
%%

coeff = 4;
level = 3;
dim = 2^level;
st = (dim+1)/2;
%%
% Different versions of this image:
%%
% 1) Just level 3 information
%
%create figure window appropriately sized
f1 = figure('Position', [100,100, 1040, 400], 'WindowStyle', 'normal', 'Color', [1, 1, 1]);

%create axes to show the profile of the level phases
a1 = axes('Units', 'pixels', 'position', [690, 50, 300, 300]); hold on; axis([5 49 -pi pi]);

step_offset1 = 1;

%plot level 3 phase
plot(1:23, angle(dt{level}(coeff,1:23)), 'r', 'LineWidth', 2.0);
plot(24:35, angle(dt{level}(coeff,24:35)), 'r', 'LineWidth', 2.0);
plot(36:64, angle(dt{level}(coeff,36:64)), 'r', 'LineWidth', 2.0);
level3_phase = angle(dt{level}(coeff,1));
theta3 = text(step_offset1, level3_phase, '\theta_3');

title('DT-CWT coefficient phase (\theta) vs step offset (x)');
xlabel('Offset of step from coefficient  (x)')
ylabel('DT-CWT phase (\theta)');

set(a1, 'Xtick', st+(0:7)*dim, 'Xticklabel', (-3:4)*dim, 'Ytick', [-pi 0 pi], 'Yticklabel', {'-Pi', '0', 'Pi'});

%Create axes to show phase arrow rotating
a2 = axes('Units', 'pixels', 'position', [360, 50, 300, 300]); hold on; axis([-1 1 -1 1]); axis off;

%Phase arrow for level 3
level3_cos = cos(level3_phase);
level3_sin = sin(level3_phase);
quiver(0, 0, level3_cos, level3_sin, 'r', 'LineWidth', 2.0,...
    'UDataSource', 'level3_cos',...
    'VDataSource', 'level3_sin');

legend('Level 3 phase \theta_3', 'position', [340 13 154 56]);
title('Complex phase of DT-CWT coefficient');

%Create axes for moving step
a3 = axes('Units', 'pixels', 'position', [50, 50, 300, 300]); hold on; axis([5 49 0.8 2.2]);
x3 = 1:64;
y3 = step_edges1D(1:64,1);
plot(x3, y3, 'k', 'LineWidth', 2.0,...
    'XDataSource', 'x3',...
    'YDataSource', 'y3');

plot([st+2*dim; st+2*dim], [0.8; 2.2], 'b:');
plot([st+4*dim; st+4*dim], [0.8; 2.2], 'b:');
plot([st+3*dim; st+3*dim], [0.8; 2.2], 'r:');

title('Positive 1D step translated past a level 3 DT-CWT coefficient');
xlabel('Offset of step from coefficient  (x)')
set(a3, 'Xtick', st+(0:7)*dim, 'Xticklabel', (-3:4)*dim, 'Yticklabel', []);
%
for offset = 5:49
    
    %Update the points that will move in the animation
    step_offset1 = offset;
    level3_phase = angle(dt{level}(coeff,offset));
    level3_cos = cos(level3_phase);
    level3_sin = sin(level3_phase);

    
    y3 = step_edges1D(1:64,offset);
    pause(.1);
    set(theta3, 'position', [step_offset1 level3_phase]);
    refreshdata(a1, 'caller');
    refreshdata(a2, 'caller');
    refreshdata(a3, 'caller');
    
    frame1 = getframe(f1);
    gif1 = frame2im(frame1);
    [gif1a map] = rgb2ind(gif1, 2^16);
    
    if offset == 5
        imwrite(gif1a, map, 'K:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\phase_3.gif', 'gif', 'WriteMode', 'overwrite', 'Loop', Inf, 'DelayTime', 0.5);
    else
        imwrite(gif1a, map, 'K:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\phase_3.gif', 'gif', 'WriteMode', 'append', 'DelayTime', 0.5);
    end
    
end
%close(f1);
%%
% 2) Level 3 and 4 information
%
%create figure window appropriately sized
f1 = figure('Position', [100,100, 1040, 400], 'WindowStyle', 'normal', 'Color', [1, 1, 1]);

%create axes to show the profile of the level phases
a1 = axes('Units', 'pixels', 'position', [690, 50, 300, 300]); hold on; axis([5 49 -pi pi]);

step_offset1 = 1;

%plot level 4 phase
plot(5:19, angle(dt{level+1}(coeff/2,1:15)), 'g:', 'LineWidth', 2.0);
plot(20:42, angle(dt{level+1}(coeff/2,16:38)), 'g:', 'LineWidth', 2.0);
plot(43:64, angle(dt{level+1}(coeff/2,39:60)), 'g:', 'LineWidth', 2.0);
level4_phase = angle(dt{level+1}(coeff/2,5));
theta4 = text(step_offset1, level4_phase, '\theta_4');

%plot level 3 phase
plot(1:23, angle(dt{level}(coeff,1:23)), 'r', 'LineWidth', 2.0);
plot(24:35, angle(dt{level}(coeff,24:35)), 'r', 'LineWidth', 2.0);
plot(36:64, angle(dt{level}(coeff,36:64)), 'r', 'LineWidth', 2.0);
level3_phase = angle(dt{level}(coeff,1));
theta3 = text(step_offset1, level3_phase, '\theta_3');

title('DT-CWT coefficient phase (\theta) vs step offset (x)');
xlabel('Offset of step from coefficient  (x)')
ylabel('DT-CWT phase (\theta)');

set(a1, 'Xtick', st+(0:7)*dim, 'Xticklabel', (-3:4)*dim, 'Ytick', [-pi 0 pi], 'Yticklabel', {'-Pi', '0', 'Pi'});

%Create axes to show phase arrow rotating
a2 = axes('Units', 'pixels', 'position', [360, 50, 300, 300]); hold on; axis([-1 1 -1 1]); axis off;

%Phase arrow for level 3
level3_cos = cos(level3_phase);
level3_sin = sin(level3_phase);
quiver(0, 0, level3_cos, level3_sin, 'r', 'LineWidth', 2.0,...
    'UDataSource', 'level3_cos',...
    'VDataSource', 'level3_sin');

%Phase arrow for level 4
level4_cos = cos(level4_phase);
level4_sin = sin(level4_phase);
quiver(0, 0, level4_cos, level4_sin, 'g:', 'LineWidth', 2.0,...
    'UDataSource', 'level4_cos',...
    'VDataSource', 'level4_sin');

legend('Level 3 phase \theta_3', 'Level 4 phase \theta_4', 'position', [340 13 154 56]);
title('Complex phase of DT-CWT coefficient');

%Create axes for moving step
a3 = axes('Units', 'pixels', 'position', [50, 50, 300, 300]); hold on; axis([5 49 0.8 2.2]);
x3 = 1:64;
y3 = step_edges1D(1:64,1);
plot(x3, y3, 'k', 'LineWidth', 2.0,...
    'XDataSource', 'x3',...
    'YDataSource', 'y3');

plot([st+2*dim; st+2*dim], [1.1; 1.9], 'b:');
plot([st+4*dim; st+4*dim], [1.1; 1.9], 'b:');
plot([st+1*dim; st+1*dim], [0.8; 2.2], 'b:');
plot([st+5*dim; st+5*dim], [0.8; 2.2], 'b:');
plot([st+3*dim; st+3*dim], [0.8; 2.2], 'r:');

title('Positive 1D step translated past a level 3 DT-CWT coefficient');
xlabel('Offset of step from coefficient  (x)')
set(a3, 'Xtick', st+(0:7)*dim, 'Xticklabel', (-3:4)*dim, 'Yticklabel', []);
%
for offset = 5:49
    
    %Update the points that will move in the animation
    step_offset1 = offset;
    level3_phase = angle(dt{level}(coeff,offset));
    level4_phase = angle(dt{level+1}(coeff/2,offset-4));
    level4_phase2 = angle(dt{level+1}(coeff/2,offset-4).^2);
    level3_cos = cos(level3_phase);
    level3_sin = sin(level3_phase);
    level4_cos = cos(level4_phase);
    level4_sin = sin(level4_phase);
    
    y3 = step_edges1D(1:64,offset);
    pause(.1);
    set(theta3, 'position', [step_offset1 level3_phase]);
    set(theta4, 'position', [step_offset1 level4_phase]);
    
    refreshdata(a1, 'caller');
    refreshdata(a2, 'caller');
    refreshdata(a3, 'caller');
    
    frame1 = getframe(f1);
    gif1 = frame2im(frame1);
    [gif1a map] = rgb2ind(gif1, 2^16);
    
    if offset == 5
        imwrite(gif1a, map, 'K:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\phase_4.gif', 'gif', 'WriteMode', 'overwrite', 'Loop', Inf, 'DelayTime', 0.5);
    else
        imwrite(gif1a, map, 'K:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\phase_4.gif', 'gif', 'WriteMode', 'append', 'DelayTime', 0.5);
    end
    
end
%close(f1);
%%
% 3) Level 3, 4 and phase-doubled information
%
%create figure window appropriately sized
f1 = figure('Position', [100,100, 1040, 400], 'WindowStyle', 'normal', 'Color', [1, 1, 1]);

%create axes to show the profile of the level phases
a1 = axes('Units', 'pixels', 'position', [690, 50, 300, 300]); hold on; axis([5 49 -pi pi]);

step_offset1 = 1;

%plot level 4 phase
% plot(5:64, angle(dt{level+1}(coeff/2,1:60)), 'g:', 'LineWidth', 2.0);
level4_phase = angle(dt{level+1}(coeff/2,5));
% theta4 = text(step_offset1, level4_phase, '\theta_4');

%plot level 3 phase
plot(1:23, angle(dt{level}(coeff,1:23)), 'r', 'LineWidth', 2.0);
plot(24:35, angle(dt{level}(coeff,24:35)), 'r', 'LineWidth', 2.0);
plot(36:64, angle(dt{level}(coeff,36:64)), 'r', 'LineWidth', 2.0);
level3_phase = angle(dt{level}(coeff,1));
theta3 = text(step_offset1, level3_phase, '\theta_3');

%plot level 4 phase doubled
level4_phase2 = angle(dt{level+1}(coeff/2,5).^2);
theta4_2 = text(step_offset1, level4_phase, '2\theta_4');
plot(5:11, angle(dt{level+1}(coeff/2,1:7).^2), 'g', 'LineWidth', 2.0);
plot(12:25, angle(dt{level+1}(coeff/2,8:21).^2), 'g', 'LineWidth', 2.0);
plot(26:36, angle(dt{level+1}(coeff/2,22:32).^2), 'g', 'LineWidth', 2.0);
plot(37:64, angle(dt{level+1}(coeff/2,33:60).^2), 'g', 'LineWidth', 2.0);

title('DT-CWT coefficient phase (\theta) vs step offset (x)');
xlabel('Offset of step from coefficient  (x)')
ylabel('DT-CWT phase (\theta)');

set(a1, 'Xtick', st+(0:7)*dim, 'Xticklabel', (-3:4)*dim, 'Ytick', [-pi 0 pi], 'Yticklabel', {'-Pi', '0', 'Pi'});

%Create axes to show phase arrow rotating
a2 = axes('Units', 'pixels', 'position', [360, 50, 300, 300]); hold on; axis([-1 1 -1 1]); axis off;

%Phase arrow for level 3
level3_cos = cos(level3_phase);
level3_sin = sin(level3_phase);
quiver(0, 0, level3_cos, level3_sin, 'r', 'LineWidth', 2.0,...
    'UDataSource', 'level3_cos',...
    'VDataSource', 'level3_sin');

%Phase arrow for level 4
level4_cos = cos(level4_phase);
level4_sin = sin(level4_phase);
quiver(0, 0, level4_cos, level4_sin, 'g:', 'LineWidth', 2.0,...
    'UDataSource', 'level4_cos',...
    'VDataSource', 'level4_sin');

%Phase-doubled arrow for level 4
level4_cos2 = cos(level4_phase2);
level4_sin2 = sin(level4_phase2);
quiver(0, 0, level4_cos2, level4_sin2, 'g', 'LineWidth', 2.0,...
    'UDataSource', 'level4_cos2',...
    'VDataSource', 'level4_sin2');

legend('Level 3 phase \theta_3', 'Level 4 phase \theta_4', 'Level 4 phase-doubled 2\theta_4', 'position', [340 13 154 56]);
title('Complex phase of DT-CWT coefficient');

%Create axes for moving step
a3 = axes('Units', 'pixels', 'position', [50, 50, 300, 300]); hold on; axis([5 49 0.8 2.2]);
x3 = 1:64;
y3 = step_edges1D(1:64,1);
plot(x3, y3, 'k', 'LineWidth', 2.0,...
    'XDataSource', 'x3',...
    'YDataSource', 'y3');

plot([st+2*dim; st+2*dim], [1.1; 1.9], 'b:');
plot([st+4*dim; st+4*dim], [1.1; 1.9], 'b:');
plot([st+1*dim; st+1*dim], [0.8; 2.2], 'b:');
plot([st+5*dim; st+5*dim], [0.8; 2.2], 'b:');
plot([st+3*dim; st+3*dim], [0.8; 2.2], 'r:');

title('Positive 1D step translated past a level 3 DT-CWT coefficient');
xlabel('Offset of step from coefficient  (x)')
set(a3, 'Xtick', st+(0:7)*dim, 'Xticklabel', (-3:4)*dim, 'Yticklabel', []);
%
for offset = 5:49
    
    %Update the points that will move in the animation
    step_offset1 = offset;
    level3_phase = angle(dt{level}(coeff,offset));
    level4_phase = angle(dt{level+1}(coeff/2,offset-4));
    level4_phase2 = angle(dt{level+1}(coeff/2,offset-4).^2);
    level3_cos = cos(level3_phase);
    level3_sin = sin(level3_phase);
    level4_cos = cos(level4_phase);
    level4_sin = sin(level4_phase);
    level4_cos2 = cos(level4_phase2);
    level4_sin2 = sin(level4_phase2);
    
    y3 = step_edges1D(1:64,offset);
    pause(.1);
    set(theta3, 'position', [step_offset1 level3_phase]);
    %set(theta4, 'position', [step_offset1 level4_phase]);
    set(theta4_2, 'position', [step_offset1 level4_phase2]);
    refreshdata(a1, 'caller');
    refreshdata(a2, 'caller');
    refreshdata(a3, 'caller');
    
    frame1 = getframe(f1);
    gif1 = frame2im(frame1);
    [gif1a map] = rgb2ind(gif1, 2^16);
    
    if offset == 5
        imwrite(gif1a, map, 'K:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\phase_42.gif', 'gif', 'WriteMode', 'overwrite', 'Loop', Inf, 'DelayTime', 0.5);
    else
        imwrite(gif1a, map, 'K:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\phase_42.gif', 'gif', 'WriteMode', 'append', 'DelayTime', 0.5);
    end
    if offset == 32
        imwrite(gif1a, map, 'K:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\phase_42_fixed.gif', 'gif', 'WriteMode', 'overwrite', 'Loop', Inf, 'DelayTime', 0.5);
    end
end
%close(f1);
%%
