%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%------- Script to build figures for Monday meeting presentation
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%%
%--------------------------------------------------------------------------
% ILP theory figures
%--------------------------------------------------------------------------

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
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% 2) ILP on block images
%--------------------------------------------------------------------------
pos_edge_128 = [100*ones(128,64) 200*ones(128,64)];
imwrite(uint8(pos_edge_128), ['K:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\pos_edge.bmp']);
dt_pos_edge = dtwavexfm2(pos_edge_128, 5, 'near_sym_b','qshift_b');
[ilp_pos_edge icp_pos_edge] = mb_dual_tree_transform(dt_pos_edge);

pos_line_128 = 100*ones(128,128); pos_line_128(:,64:65) = 200;
imwrite(uint8(pos_line_128), ['K:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\pos_line.bmp']);
dt_pos_line = dtwavexfm2(pos_line_128, 5, 'near_sym_b','qshift_b');
[ilp_pos_line icp_pos_line] = mb_dual_tree_transform(dt_pos_line);

neg_line_128 = 200*ones(128,128); neg_line_128(:,64:65) = 100;
imwrite(uint8(neg_line_128), ['K:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\neg_line.bmp']);
dt_neg_line = dtwavexfm2(neg_line_128, 5, 'near_sym_b','qshift_b');
[ilp_neg_line icp_neg_line] = mb_dual_tree_transform(dt_neg_line);

f1 = figure;
for level = 2:4
    dim = 2^level;
    st = (dim+1)/2;

    %[max_icp_pos] = max(icp_pos_edge{level}, [], 3);
    
    
    %ILP of pos edge
    [max_ilp_pos_edge] = max(ilp_pos_edge{level}, [], 3);
    f = figure('Position', [100,100, 512, 512], 'WindowStyle', 'normal', 'Color', [1, 1, 1]);
    axes('Units', 'pixels', 'position', [0, 0, 512, 512]); axis image; hold on;
    image(pos_edge_128); axis image; colormap(gray(256));
    quiver(st:dim:129-st, st:dim:129-st, real(max_ilp_pos_edge), -imag(max_ilp_pos_edge), 'r');
    saveas(f, ['K:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\ilp_edge', num2str(level), '.png']);
    close(f);
    imwrite(complex2rgb(max_ilp_pos_edge), ['K:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\ilp_edge_mp', num2str(level), '.png']);
    
    %ILP of pos line
    [max_ilp_pos_line] = max(ilp_pos_line{level}, [], 3);
    f = figure('Position', [100,100, 512, 512], 'WindowStyle', 'normal', 'Color', [1, 1, 1]);
    axes('Units', 'pixels', 'position', [0, 0, 512, 512]); axis image; hold on;
    image(pos_line_128); axis image; colormap(gray(256)); hold on;
    quiver(st:dim:129-st, st:dim:129-st, real(max_ilp_pos_line), -imag(max_ilp_pos_line), 'r');
    saveas(f, ['K:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\ilp_pos_line', num2str(level), '.png']);
    close(f);
    imwrite(complex2rgb(max_ilp_pos_line), ['K:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\ilp_pos_line_mp', num2str(level), '.png']);
    
    %ILP of neg line
    [max_ilp_neg_line] = max(ilp_neg_line{level}, [], 3);
    f = figure('position', [100,100, 512, 512], 'WindowStyle', 'normal', 'Color', [1, 1, 1]);
    axes('Units', 'pixels', 'position', [0, 0, 512, 512]); axis image; hold on;
    image(neg_line_128); axis image; colormap(gray(256));
    quiver(st:dim:129-st, st:dim:129-st, real(max_ilp_neg_line), -imag(max_ilp_neg_line), 'r');
    saveas(f, ['K:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\ilp_neg_line', num2str(level), '.png']);
    close(f);
    imwrite(complex2rgb(max_ilp_neg_line), ['K:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\ilp_neg_line_mp', num2str(level), '.png']);
    
    %ICP of pos edge
    [max_icp_pos_edge] = max(icp_pos_edge{level}, [], 3);
    imwrite(complex2rgb(abs(max_icp_pos_edge).*exp(i*2*angle(max_icp_pos_edge))), ['K:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\icp_edge_mp', num2str(level), '.png']);
    
    %ICP of pos line
    [max_icp_pos_line] = max(icp_pos_line{level}, [], 3);
    imwrite(complex2rgb(abs(max_icp_pos_line).*exp(i*2*angle(max_icp_pos_line))), ['K:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\icp_pos_line_mp', num2str(level), '.png']);
    
    %ICP of neg line
    [max_icp_neg_line] = max(icp_neg_line{level}, [], 3);
    imwrite(complex2rgb(abs(max_icp_neg_line).*exp(i*2*angle(max_icp_neg_line))), ['K:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\icp_neg_line_mp', num2str(level), '.png']);
    
    %ILP compass plots
    figure(f2);
    subplot(1,3,1);
    compass(max_ilp_pos_line(2^(6-level), 2^(6-level)), colors(level)); hold on;
    subplot(1,3,2);
    compass(max_ilp_pos_edge(2^(6-level), 2^(6-level)), colors(level)); hold on;
    subplot(1,3,3);
    compass(max_ilp_pos_line(2^(6-level), 2^(6-level)), colors(level)); hold on;
end
%%
% Write images for all sub-bands of line
mkdir('K:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\pos_line');
pos_line = 100*ones(256); pos_line(:,128) = 200;
pos_line = imrotate(pos_line, 15, 'crop');
pos_line = pos_line(65:192, 65:192);
imwrite(uint8(pos_line), 'K:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\pos_line_12.bmp');
dt_pos_line = dtwavexfm2(pos_line, 5, 'near_sym_b','qshift_b');
[ilp_pos_line icp_pos_line] = mb_dual_tree_transform(dt_pos_line);
%%
for level = 1:4
    m_dt = max(abs(dt_pos_line{level}(:)));
    m_ilp = max(abs(ilp_pos_line{level}(:)));
    m_icp = max(abs(icp_pos_line{level}(:)));
    for band = 1:6
        imwrite(imresize(complex2rgb(dt_pos_line{level}(:,:,band), [-pi pi], m_dt), 2^level, 'nearest'),...
            ['K:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\pos_line\dt_', num2str(level), '_', num2str(band) '.bmp']);
        imwrite(imresize(complex2rgb(ilp_pos_line{level}(:,:,band), [-pi pi], m_ilp), 2^level, 'nearest'),...
            ['K:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\pos_line\ilp_', num2str(level), '_', num2str(band) '.bmp']);
        imwrite(imresize(complex2rgb(icp_pos_line{level}(:,:,band), [-pi pi], m_icp), 2^level, 'nearest'),...
            ['K:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\pos_line\icp_', num2str(level), '_', num2str(band) '.bmp']);
    end
    imwrite(imresize(complex2rgb(max(dt_pos_line{level}, [], 3), [-pi pi], m_dt), 2^level, 'nearest'),...
        ['K:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\pos_line\dt_', num2str(level), '_max.bmp']);
    imwrite(imresize(complex2rgb(max(ilp_pos_line{level}, [], 3), [-pi pi], m_ilp), 2^level, 'nearest'),...
        ['K:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\pos_line\ilp_', num2str(level), '_max.bmp']);
    imwrite(imresize(complex2rgb(max(icp_pos_line{level}, [], 3), [-pi pi], m_icp), 2^level, 'nearest'),...
        ['K:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\pos_line\icp_', num2str(level), '_max.bmp']);
end
%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Create radial lines image
%--------------------------------------------------------------------------
angles = 0:18:162;
pos_line512 = zeros(512); pos_line512([17:224 289:496], 256:257) = 1;
radial_lines = zeros(512);
for ii = 1:length(angles)
    radial_lines = radial_lines + imrotate(pos_line512, angles(ii), 'crop');
end
radial_lines = 100*(radial_lines+1);
figure; image(radial_lines); axis image; colormap(gray(256));
%imwrite(uint8(radial_lines), 'K:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\radial_lines.bmp');

dt_radial_lines = dtwavexfm2(radial_lines, 5, 'near_sym_b','qshift_b');
[ilp_radial_lines icp_radial_lines] = mb_dual_tree_transform(dt_radial_lines);
%%
for level = 1:4
    %[max_icp_radial_lines] = max(icp_radial_lines{level}, [], 3);
    [max_icp_radial_lines] = compute_max_icp(icp_radial_lines{level});
    [max_ilp_radial_lines] = max(ilp_radial_lines{level}, [], 3);
    
    figure; image(complex2rgb(abs(max_icp_radial_lines).*exp(i*2*angle(max_icp_radial_lines)))); axis image;
    figure; image(complex2rgb(max_ilp_radial_lines)); axis image;
%     imwrite(complex2rgb(max_ilp_radial_lines), ['K:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\ilp_radial_mp', num2str(level), '.bmp']);
%     imwrite(complex2rgb(abs(max_icp_radial_lines).*exp(i*2*angle(max_icp_radial_lines))), ['K:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\icp_radial_mp', num2str(level), '.bmp']);
end

%%
%--------------------------------------------------------------------------
% Create different line profiles image

%Create gaussian line profile
x = -64:1:63;
halfwidth = 3;
sigma2 = (halfwidth^2) / log(2);
ymax = 1/sqrt(2*pi*sigma2);

scaling = 100 / ymax;
g_x = scaling*exp(-(x.^2 / sigma2)) / sqrt(2*pi*sigma2);
figure; plot(x, g_x);
hold on;

%imwrite(uint8(pos_line_gauss), 'K:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\gauss_line.bmp');
%imwrite(uint8(pos_line_wide), 'K:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\wide_line.bmp');
%
dt_line = dtwavexfm2(pos_line_wide, 5, 'near_sym_b','qshift_b');
dt_gauss = dtwavexfm2(pos_line_gauss, 5, 'near_sym_b','qshift_b');

plot(0, ymax, 'gx');
plot([halfwidth halfwidth], [0 scaling*ymax], 'r');
plot([-halfwidth -halfwidth], [0 scaling*ymax], 'r');
plot([halfwidth -halfwidth], scaling*[ymax ymax]/2, 'gx');

line_gauss = repmat(g_x, 128, 1)+100;
figure; image(line_gauss); axis image; colormap(gray(256));

line_wide = 100*ones(128,128); line_wide(:,62:67) = 200;
figure; image(line_wide); axis image; colormap(gray(256));
%%
imwrite(uint8(line_gauss), 'K:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\line_gauss.bmp');
imwrite(uint8(line_wide), 'K:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\line_wide.bmp');
%%
dt_line = dtwavexfm2(line_wide, 5, 'near_sym_b','qshift_b');
dt_gauss = dtwavexfm2(line_gauss, 5, 'near_sym_b','qshift_b');

[ilp_line icp_line] = mb_dual_tree_transform(dt_line);
[ilp_gauss icp_gauss] = mb_dual_tree_transform(dt_gauss);
%%
f1 = figure; 
subplot(1,2,1);
image(line_wide); axis image; colormap(gray(256)); hold on;
subplot(1,2,2);
image(line_gauss); axis image; colormap(gray(256)); hold on;

f2 = figure; 

colors = 'rgbmy';
r = 62; c = 64;
for level = 1:5
    
    dim = 2^level;
    st = (dim+1)/2;
    rl = ceil(r / 2^level);
    cl = ceil(c / 2^level);
    
    [max_ilp_line] = max(ilp_line{level}, [], 3);
    [max_ilp_gauss] = max(ilp_gauss{level}, [], 3);
    
    figure(f1);
    subplot(1,2,1);        
    quiver(st:dim:129-st, st:dim:129-st, real(max_ilp_line), -imag(max_ilp_line), colors(level));
    subplot(1,2,2);
    quiver(st:dim:129-st, st:dim:129-st, real(max_ilp_gauss), -imag(max_ilp_gauss), colors(level));
    
    figure(f2); hold on;
    subplot(1,2,1); axis equal; hold on; axis([-100 5000 -100 5000])
    compass(max_ilp_line(rl, cl), colors(level));
    subplot(1,2,2); axis equal; hold on; axis([-100 5000 -100 5000])
    compass(max_ilp_gauss(rl, cl), [colors(level), ':']);
end
%%
% Display line profiles
f1 = figure('Position', [100,100, 192, 128], 'WindowStyle', 'normal', 'Color', [1, 1, 1]);
a1 = axes('Units', 'pixels', 'position', [0, 0, 192, 128]); hold on;
plot(1:128, line_wide(1,:), 'k', 'linewidth', 2.0);  axis([1 128 50 250]); axis off;
saveas(f1, 'K:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\line_wide_profile.bmp');
f2 = figure('Position', [100,100, 192, 128], 'WindowStyle', 'normal', 'Color', [1, 1, 1]);
a2 = axes('Units', 'pixels', 'position', [0, 0, 192, 128]); hold on;
plot(1:128, line_gauss(1,:), 'k', 'linewidth', 2.0); axis([1 128 50 250]); axis off;
saveas(f2, 'K:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\line_gauss_profile.bmp');
%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Create ILP/ICP images for real mammogram plots, can we show a normal
% patch and a radial patch?
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%Find streaky patch from ILP report

%try backgrounds, 006, 008 and 004 for radial texture
bg = double(imread('C:\isbe\dev\background\images\mass\mass008.bmp'));
bg = bg(142:1100, :);
dt_bg = dtwavexfm2(bg, 7);

[ilp_bg icp_bg] = mb_dual_tree_transform(dt_bg);
for level = 2:5
    
    [max_ilp{level}] = max(ilp_bg{level}, [], 3);
    [max_icp{level}] = max(icp_bg{level}, [], 3);
    
    ilp_mag = histeq(abs(max_ilp{level}) / max(abs(max_ilp{level}(:))), 1 ./ (1:256));
    icp_mag = histeq(abs(max_icp{level}) / max(abs(max_icp{level}(:))), 1 ./ (1:256));
    temp = icp_mag.*exp(i*angle(max_icp{level}));
    figure; title(['ICP level ', num2str(level)]);
    image(complex2rgb(ilp_mag.*exp(i*2*angle(max_icp{level})))); axis image; hold on;
    %quiver(real(max_icp{level}), -imag(max_icp{level}), 'w');
    %figure;
    %quiver(real(temp), -imag(temp)); axis image;
    figure; title(['ILP level ', num2str(level)]);
    image(complex2rgb(ilp_mag.*exp(i*angle(max_ilp{level})))); axis image;
    hold on;
   % quiver(real(max_ilp{level}), -imag(max_ilp{level}), 'w');
   
    imwrite(complex2rgb(ilp_mag.*exp(i*angle(max_ilp{level}))),...
        ['K:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\ilp_bg_mp', num2str(level), '.bmp']);
    imwrite(complex2rgb(ilp_mag.*exp(i*2*angle(max_icp{level}))),...
        ['K:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\icp_bg_mp', num2str(level), '.bmp']);
end
%%

write_im_from_colormap(bg, 'C:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\mass008.bmp', gray(256));
%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Streaky background
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
streaky = rot90(double(imread('C:\isbe\dev\background\images\normal1024\o04_011RML_1024_3309_2896.bmp')));
streaky_dt = dtwavexfm2(streaky, 7);

[streaky_ilp streaky_icp] = mb_dual_tree_transform(streaky_dt);
%%
for level = 2:5
    
    [max_ilp{level}] = max(streaky_ilp{level}, [], 3);
    [max_icp{level}] = max(streaky_icp{level}, [], 3);
    
    ilp_mag = histeq(abs(max_ilp{level}) / max(abs(max_ilp{level}(:))), 1 ./ (1:256));
    icp_mag = histeq(abs(max_icp{level}) / max(abs(max_icp{level}(:))), 1 ./ (1:256));
    temp = icp_mag.*exp(i*angle(max_icp{level}));
    figure; title(['ICP level ', num2str(level)]);
    image(complex2rgb(ilp_mag.*exp(i*2*angle(max_icp{level})))); axis image; hold on;
    %quiver(real(max_icp{level}), -imag(max_icp{level}), 'w');
    %figure;
    %quiver(real(temp), -imag(temp)); axis image;
    figure; title(['ILP level ', num2str(level)]);
    image(complex2rgb(ilp_mag.*exp(i*angle(max_ilp{level})))); axis image;
    hold on;
   % quiver(real(max_ilp{level}), -imag(max_ilp{level}), 'w');
   
    imwrite(complex2rgb(ilp_mag.*exp(i*angle(max_ilp{level}))),...
        ['K:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\ilp_streaky_mp', num2str(level), '.bmp']);
    imwrite(complex2rgb(ilp_mag.*exp(i*2*angle(max_icp{level}))),...
        ['K:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\icp_streaky_mp', num2str(level), '.bmp']);
end
%%
[response, orientations] = line_operator_conv(streaky, 24, 5, 15);
%%
write_im_from_colormap(streaky, 'C:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\streaky.bmp', gray(256));
%%
pectoral = double(imread('D:\isbe\dev\background\images\pectoral\pectoral010.bmp'));
pectoral_dt = dtwavexfm2(pectoral, 7);

[pectoral_ilp pectoral_icp] = mb_dual_tree_transform(pectoral_dt);
%
for level = 2:5
    
    [max_ilp{level}] = max(pectoral_ilp{level}, [], 3);
    [max_icp{level}] = max(pectoral_icp{level}, [], 3);
    
    ilp_mag = histeq(abs(max_ilp{level}) / max(abs(max_ilp{level}(:))), 1 ./ (1:256));
    icp_mag = histeq(abs(max_icp{level}) / max(abs(max_icp{level}(:))), 1 ./ (1:256));
    temp = icp_mag.*exp(i*angle(max_icp{level}));
    figure; title(['ICP level ', num2str(level)]);
    image(complex2rgb(ilp_mag.*exp(i*2*angle(max_icp{level})))); axis image; hold on;
    %quiver(real(max_icp{level}), -imag(max_icp{level}), 'w');
    %figure;
    %quiver(real(temp), -imag(temp)); axis image;
    figure; title(['ILP level ', num2str(level)]);
    image(complex2rgb(ilp_mag.*exp(i*angle(max_ilp{level})))); axis image;
    hold on;
   % quiver(real(max_ilp{level}), -imag(max_ilp{level}), 'w');
   
    imwrite(complex2rgb(ilp_mag.*exp(i*angle(max_ilp{level}))),...
        ['K:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\ilp_pectoral_mp', num2str(level), '.bmp']);
    imwrite(complex2rgb(ilp_mag.*exp(i*2*angle(max_icp{level}))),...
        ['K:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\icp_pectoral_mp', num2str(level), '.bmp']);
end

write_im_from_colormap(pectoral, 'C:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\pectoral.bmp', gray(256));
%%
% Pectoral image and synthetic pectoral image
%pec = double((imread('C:\isbe\dev\background\images\pectoral\pectoral07.bmp')));
pec = double(rgb2gray(imread('C:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\synthetic_pec.bmp')));
figure; imagesc(pec); axis image; colormap(gray(256));
%
dt_pec = dtwavexfm2(pec, 6);
[pec_ilp pec_icp] = mb_dual_tree_transform(dt_pec);

% w = [-3 -1; -3 -3; -1 -3; 1 -3; 3 -3; 3 -1]*pi/2.15; % Nominally j * pi/2, but reduced a bit due to asymmetry of subband freq responses.
% p = [1 3]/4;  % Interpolation points
%
% ilp_all = pec_ilp{5};
% icp_all = pec_icp{5};
%    
for level = 4:-1:1
    
    [max_ilp{level}] = max(pec_ilp{level}, [], 3);
    [max_icp{level}] = max(pec_icp{level}, [], 3);
%     [r c] = size(pec_ilp{level}(:,:,1));
%     for band = 1:6,
%         % clipping extra space, for the moment
%         %upscale_ilp = cpxinterp2(ilp_all(:,:,band), p-0.5, w(band,:),'spline');
%         upscale_ilp = kron(ilp_all(:,:,band), ones(2));
%         temp_ilp(:,:,band) = max(upscale_ilp(1:r, 1:c), pec_ilp{level}(:,:,band));
%         
%         %upscale_icp = cpxinterp2(icp_all(:,:,band), p-0.5, w(band,:),'spline'); 
%         upscale_icp = kron(icp_all(:,:,band), ones(2));
%         temp_icp(:,:,band) = max(upscale_icp(1:r, 1:c), pec_icp{level}(:,:,band));
%     end
%     ilp_all = temp_ilp;
%     icp_all = temp_icp;
%     clear temp*
    ilp_mag = histeq(abs(max_ilp{level}) / max(abs(max_ilp{level}(:))), 1 ./ (1:256));
    icp_mag = histeq(abs(max_icp{level}) / max(abs(max_icp{level}(:))), 1 ./ (1:256));
    temp = icp_mag.*exp(i*angle(max_icp{level}));
    figure;
    image(complex2rgb(icp_mag.*exp(i*2*angle(max_icp{level})))); axis image;
    %figure;
    %quiver(real(temp), -imag(temp)); axis image;
    figure;
    image(complex2rgb(ilp_mag.*exp(i*angle(max_ilp{level})))); axis image;
    
%     imwrite(complex2rgb(ilp_mag.*exp(i*angle(max_ilp{level}))),...
%         ['C:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\ilp_streaky_mp', num2str(level), '.bmp']);
%     imwrite(complex2rgb(icp_mag.*exp(i*angle(max_icp{level}))),...
%         ['C:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\icp_streaky_mp', num2str(level), '.bmp']);
end
% [max_ilp_all] = max(ilp_all, [], 3);
% [max_icp_all] = max(icp_all, [], 3);
% figure; image(complex2rgb(max_ilp_all)); axis image;
% figure; image(complex2rgb(abs(max_icp_all).*exp(i*2*angle(max_icp_all)))); axis image;
%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%-Tree image
tree = (double(rgb2gray(imread('K:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\tree3.JPG'))));
imwrite(uint8(tree),...
        'K:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\tree\tree_gray.bmp');
figure; imagesc(tree); axis image; colormap(gray(256));
[dt_tree sc_tree] = dtwavexfm2(tree, 6);
[ilp_tree icp_tree] = mb_dual_tree_transform(dt_tree); 
%%
mkdir('K:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\tree');
for level = 1:5
    
    [max_ilp] = max(ilp_tree{level}, [], 3);
    [max_icp] = max(icp_tree{level}, [], 3);
    [max_dt] = max(dt_tree{level}, [], 3);
    
    ilp_mag = histeq(abs(max_ilp) / max(abs(max_ilp(:))), 1 ./ (1:256));
    icp_mag = histeq(abs(max_icp) / max(abs(max_icp(:))), 1 ./ (1:256));
    figure('name', ['ICP level ', num2str(level)]);
    image(complex2rgb(icp_mag.*exp(i*2*angle(max_icp)))); axis image;
    imwrite(complex2rgb(icp_mag.*exp(i*2*angle(max_icp))),...
        ['K:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\tree\tree_icp_', num2str(level), '.bmp']);
    
    figure('name', ['ILP level ', num2str(level)]);
    image(complex2rgb(ilp_mag.*exp(i*angle(max_ilp)))); axis image;
    imwrite(complex2rgb(ilp_mag.*exp(i*angle(max_ilp))),...
        ['K:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\tree\tree_ilp_', num2str(level), '.bmp']);
    
    figure('name', ['DT level ', num2str(level)]);
    imagesc(histeq(abs(max_dt) / max(abs(max_dt(:))), 1 ./ (1:256))); colormap(gray(256)); axis image;
    imwrite(histeq(abs(max_dt) / max(abs(max_dt(:))), 1 ./ (1:256)),...
        ['K:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\tree\tree_dt_', num2str(level), '.bmp']);
    
    dim = 2^level;
    st = (dim+1)/2;
    [r c] = size(max_icp);
    
    figure; imagesc(tree); axis image; colormap(gray(256)); hold on;        
    quiver(st+(0:c-1)*dim, st+(0:r-1)*dim, real(icp_mag.*exp(i*angle(max_icp))), -imag(icp_mag.*exp(i*angle(max_icp))), 'm');
    
    if level == 2;
        figure;
        for band = 1:6
            abs_band = abs(dt_tree{level}(:,:,band));
            subplot(2,3,band);
            imagesc(histeq(abs_band / max(abs_band(:)), 1 ./ (1:256))); colormap(gray(256)); axis image;
            
            imwrite(histeq(abs_band / max(abs_band(:)), 1 ./ (1:256)),...
                ['K:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\tree\tree_dt_2_', num2str(band), '.bmp']);
        end
    end
    
    %icp_all_levels(:,:,level) = imresize(max_icp, size(icp_tree{1}(:,:,1)), 'nearest');
            
    
    %m_dt = max(abs(dt_tree{level}(:)));
    %m_ilp = max(abs(ilp_tree{level}(:)));
    %m_icp = max(abs(icp_tree{level}(:)));
    
    %f1 = figure;
    %f2 = figure;
    %f3 = figure;
%     for band = 1:6
%         imwrite(imresize(complex2rgb(dt_pos_line{level}(:,:,band), [-pi pi], m_dt), 2^level, 'nearest'),...
%             ['K:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\pos_line\dt_', num2str(level), '_', num2str(band) '.bmp']);
%         imwrite(imresize(complex2rgb(ilp_pos_line{level}(:,:,band), [-pi pi], m_ilp), 2^level, 'nearest'),...
%             ['K:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\pos_line\ilp_', num2str(level), '_', num2str(band) '.bmp']);
%         imwrite(imresize(complex2rgb(icp_pos_line{level}(:,:,band), [-pi pi], m_icp), 2^level, 'nearest'),...
%             ['K:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\pos_line\icp_', num2str(level), '_', num2str(band) '.bmp']);
%        figure(f1); subplot(2,3,band); imagesc(abs(dt_tree{level}(:,:,band))); colormap(gray(256)); axis image; %caxis([0 m_dt]);
%         figure(f2); subplot(2,3,band); image(complex2rgb(ilp_tree{level}(:,:,band), [-pi pi], m_ilp)); axis image
%         figure(f3); subplot(2,3,band); image(complex2rgb(icp_tree{level}(:,:,band), [-pi pi], m_icp)); axis image
%     end
end
%%
max_icp_all_levels = max(icp_all_levels(1:2:end,1:2:end,1:4), [], 3);
icp_mag = histeq(abs(max_icp_all_levels) / max(abs(max_icp_all_levels(:))), 1 ./ (1:256));
figure('name', 'ICP all levels');
image(complex2rgb(icp_mag.*exp(i*2*angle(max_icp_all_levels)))); axis image;

figure; imagesc(tree); axis image; colormap(gray(256)); hold on;
[r c] = size(max_icp_all_levels);
 quiver((3:8:8*c - 1)/2, (3:8:8*r - 1)/2, real(icp_mag.*exp(i*angle(max_icp_all_levels))), -imag(icp_mag.*exp(i*angle(max_icp_all_levels))));
%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Create Colorbar for phase magnitude plots
%--------------------------------------------------------------------------
phase_colorbar = repmat([0; exp(i*((2*pi) * (1:360)' ./ 360))], 1, 20);
imwrite(complex2rgb(phase_colorbar), 'K:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\phase_colorbar.bmp');
mag_colorbar = repmat((360:-1:0)', 1, 20);
imwrite(complex2rgb(mag_colorbar), 'K:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\mag_colorbar.bmp'); axis image
%%
% Colorwheel
[xx yy] = meshgrid(linspace(-1, 1, 256), linspace(1, -1, 256));
colorwheel = complex2rgb(complex(xx, yy), [-pi pi], 1);
mask = xx.^2 + yy.^2 > 1;
colorwheel(cat(3, mask, mask, mask)) = 1;
figure; imagesc(colorwheel); axis image; colormap(gray(256));
imwrite(colorwheel, 'C:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\colorwheel.bmp');

colorwheel2 = complex2rgb(abs(complex(xx, yy)).*exp(2*i*angle(complex(xx, yy))), [-pi pi], 1);
colorwheel2(cat(3, mask, mask, mask)) = 1;
figure; imagesc(colorwheel2); axis image; colormap(gray(256));
imwrite(colorwheel2, 'C:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\colorwheel2.bmp');
%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Phase gradient for angled lines
%--------------------------------------------------------------------------
%%
line_0 = zeros(128,128); line_0(17:112,64:65) = 1;
line_12 = imrotate(line_0, 12, 'crop');
line_24 = imrotate(line_0, 24, 'crop');
dt_line_0 = dtwavexfm2(line_0, 5, 'near_sym_b','qshift_b');
dt_line_12 = dtwavexfm2(line_12, 5, 'near_sym_b','qshift_b');
dt_line_24 = dtwavexfm2(line_24, 5, 'near_sym_b','qshift_b');

w = [-3 -1; -3 -3; -1 -3; 1 -3; 3 -3; 3 -1]*pi/2.15;
p = [-3 -1 1 3]/8;

interp_0 = cpxinterp2(dt_line_0{2}(:,:,4), p, w(4,:),'linear');
interp_12 = cpxinterp2(dt_line_12{2}(:,:,4), p, w(4,:),'linear');
interp_24 = cpxinterp2(dt_line_24{2}(:,:,4), p, w(4,:),'linear');
%
cinterp_0 = complex2rgb(interp_0);
cinterp_12 = complex2rgb(interp_12);
cinterp_24 = complex2rgb(interp_24);
figure; image(cinterp_0(33:96, 33:96, :)); axis image;
figure; image(cinterp_12(33:96, 33:96, :)); axis image;
figure; image(cinterp_24(33:96, 33:96, :)); axis image;
%%
figure;
subplot(1,3,1); image(cinterp_0(49:80, 64, :));
subplot(1,3,2); image(cinterp_12(49:80, 64, :));
subplot(1,3,3); image(cinterp_24(49:80, 64, :));
%%
figure; 
subplot(1,3,1); plot(1:32, angle(interp_0(49:80, 64))); axis([0 32 -pi pi]);
subplot(1,3,2); plot(1:32, angle(interp_12(49:80, 64))); axis([0 32 -pi pi]);
subplot(1,3,3); plot(1:32, angle(interp_24(49:80, 64))); axis([0 32 -pi pi]);
%%
figure;
subplot(1,3,1); imagesc(angle(interp_0(49:80, 64))); colormap(hsv(256)); caxis([-pi pi]);
subplot(1,3,2); imagesc(angle(interp_12(49:80, 64))); colormap(hsv(256)); caxis([-pi pi]);
subplot(1,3,3); imagesc(angle(interp_24(49:80, 64))); colormap(hsv(256)); caxis([-pi pi]);
%%
imwrite(uint8(100*(1+line_0)), 'C:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\line0.bmp');
imwrite(uint8(100*(1+line_12)), 'C:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\line12.bmp');
imwrite(uint8(100*(1+line_24)), 'C:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\line24.bmp');
imwrite(cinterp_0, 'C:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\dt_line0.bmp');
imwrite(cinterp_12, 'C:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\dt_line12.bmp');
imwrite(cinterp_24, 'C:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\dt_line24.bmp');

imwrite(cinterp_0(33:96, 33:96, :), 'C:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\dt_line0z.bmp');
imwrite(cinterp_12(33:96, 33:96, :), 'C:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\dt_line12z.bmp');
imwrite(cinterp_24(33:96, 33:96, :), 'C:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\dt_line24z.bmp');
%%
write_im_from_colormap(imresize(angle(interp_0(49:80, 64)), 16, 'nearest'),...
    'C:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\dt_line0l.bmp',...
    hsv(256), [-pi pi]);
write_im_from_colormap(imresize(angle(interp_12(49:80, 64)), 16, 'nearest'),...
    'C:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\dt_line12l.bmp',...
    hsv(256), [-pi pi]);
write_im_from_colormap(imresize(angle(interp_24(49:80, 64)), 16, 'nearest'),...
    'C:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\dt_line24l.bmp',...
    hsv(256), [-pi pi]);

%%

%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% CLS detection
%--------------------------------------------------------------------------
GaborFilterArgs.Normalise = 1;
GaborFilterArgs.HighPassFilter = true;

mammogram = bg008(142:1100,:);

%build Gaussian Pyramid of mammogram
[gp p_sizes] = buildGpyr(mammogram, 4);
gp = mb_change_pyramid_form(gp, p_sizes, 'g'); clear p_sizes

%Perform CLS detection on each level of the pyramid
cls = cell(4,1);

for level = 2:4
    cls{level} = mb_cls_selection('ImageIn', gp{level},...
        'GaborFilterArgs', GaborFilterArgs,...
        'MinLength', 12, 'Connectivity', 0, 'Thin', 1,...
        'GaborThreshold', 0,...
        'GradientThreshold', 0, 'GradientAlignment', pi/12,...
        'NMS', true, 'Debug', false);
end
%%
%--------------------------------------------------------------------------
%Phase-magnitude swap - lenna/einstein
%--------------------------------------------------------------------------
einstein = imread('C:\isbe\matlab_code\trunk\matlabPyrTools\einstein.pgm');
lenna = u_load('C:\isbe\matlab_code\trunk\wavelet_dtcwt_toolbox4_1\lenna');

einstein_fft = fft2(double(einstein));
lenna_fft = fft2(double(lenna));
einstein2 = ifft2(einstein_fft);


el = real(ifft2(abs(einstein_fft).* exp(i*angle(lenna_fft))));
le = real(ifft2(exp(i*angle(einstein_fft)) .* abs(lenna_fft)));
figure; imagesc(el); axis image; colormap(gray(256));
figure; imagesc(le); axis image; colormap(gray(256));
figure; imagesc(el); axis image; colormap(gray(256));

imwrite(einstein, 'K:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\einstein.bmp'); axis image
imwrite(uint8(lenna), 'K:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\lenna.bmp'); axis image
write_im_from_colormap(el, 'K:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\el.bmp', gray(256));
write_im_from_colormap(le, 'K:\isbe\conferences_and_symposia\Monday Meeting 09-03-23\figures\le.bmp', gray(256)); 