load C:\isbe\dev\files\small_files.mat
load C:\isbe\dev\files\large_files.mat
%%
%
% Showing outlines of small and large masses
%

f_small1 = figure('name', 'Small masses: 1 to 25');
f_large1 = figure('name', 'Large masses: 1 to 25');
f_small2 = figure('name', 'Small masses: 26 to 50');
f_large2 = figure('name', 'Large masses: 26 to 50');

for ii = 1:25
    % small masses 1 to 25
    load(['C:\isbe\dev\masses\', small_files(ii).name]);
    idx = round(linspace(1, length(mass.mass_outline(:,1)), 501));
    idx(end) = [];
    x_shape = mass.mass_outline(idx,:);
    figure(f_small1); subplot(5,5,ii);
    plot(x_shape([1:end,1],1), x_shape([1:end,1],2));
    clear mass;
    
    % large masses 1 to 25
    load(['C:\isbe\dev\masses\', large_files(ii).name]);
    idx = round(linspace(1, length(mass.mass_outline(:,1)), 501));
    idx(end) = [];
    x_shape = mass.mass_outline(idx,:);
    figure(f_large1); subplot(5,5,ii);
    plot(x_shape([1:end,1],1), x_shape([1:end,1],2));
    clear mass;
    
    % small masses 26 to 50
    load(['C:\isbe\dev\masses\', small_files(ii+25).name]);
    idx = round(linspace(1, length(mass.mass_outline(:,1)), 501));
    idx(end) = [];
    x_shape = mass.mass_outline(idx,:);
    figure(f_small2); subplot(5,5,ii);
    plot(x_shape([1:end,1],1), x_shape([1:end,1],2));
    clear mass;
    
    % large masses 26 to 50
    load(['C:\isbe\dev\masses\', large_files(ii+25).name]);
    idx = round(linspace(1, length(mass.mass_outline(:,1)), 501));
    idx(end) = [];
    x_shape = mass.mass_outline(idx,:);
    figure(f_large2); subplot(5,5,ii);
    plot(x_shape([1:end,1],1), x_shape([1:end,1],2));
    clear mass;
end
%%
%
% Plotting procrustes aligned shapes, large and small
%
colours = 'rgbymck';

%small shapes
X_shape_s = zeros(50, 1000);
mass_areas_s = zeros(50,1);

for ii = 1:50
    load(['C:\isbe\dev\masses\', small_files(ii).name]);
    idx = round(linspace(1, length(mass.mass_outline(:,1)), 501));
    idx(end) = [];
    X_shape_s(ii,:) = [mass.mass_outline(idx,1)', mass.mass_outline(idx,2)'];    
    mass_areas_s(ii) = mass.mass_area;
    clear mass;
end

[mean_X, P_x, B_x, L_x,...
    mean_scale, P_scale, B_scale, L_scale mean_target, shape_scale, X_pro_s]...
    = shape_model(X_shape_s, mean(mass_areas_s), 0.98);
clear mean_X P_x B_x L_x mean_scale P_scale B_scale L_scale mean_target shape_scale;

figure('name', 'procrustes aligned small shapes'); hold on;
for ii = 1:50;
    plot(X_pro_s(ii,[1:500,1],1), X_pro_s(ii, [501:end,501]), colours(rem(ii,7)+1));
end

%large shapes
X_shape_l = zeros(50, 1000);
mass_areas_l = zeros(50,1);

for ii = 1:51
    load(['C:\isbe\dev\masses\', large_files(ii).name]);
    idx = round(linspace(1, length(mass.mass_outline(:,1)), 501));
    idx(end) = [];
    X_shape_l(ii,:) = [mass.mass_outline(idx,1)', mass.mass_outline(idx,2)'];    
    mass_areas_l(ii) = mass.mass_area;
    clear mass;
end

[mean_X, P_x, B_x, L_x,...
    mean_scale, P_scale, B_scale, L_scale mean_target, shape_scale, X_pro_l]...
    = shape_model(X_shape_l, mean(mass_areas_l), 0.98);
clear mean_X P_x B_x L_x mean_scale P_scale B_scale L_scale mean_target shape_scale;

figure('name', 'procrustes aligned large shapes'); hold on;
for ii = 1:51;
    plot(X_pro_l(ii,[1:500,1],1), X_pro_l(ii, [501:end,501]), colours(rem(ii,7)+1));
end
%%
%all shapes
X_shape = zeros(50, 1000);
mass_areas = zeros(50,1);

for ii = 1:101
    load(['C:\isbe\dev\masses\', u_files1(ii).name]);
    idx = round(linspace(1, length(mass.mass_outline(:,1)), 501));
    idx(end) = [];
    X_shape(ii,:) = [mass.mass_outline(idx,1)', mass.mass_outline(idx,2)'];    
    mass_areas(ii) = mass.mass_area;
    clear mass;
end

[mean_X, P_x, B_x, L_x,...
    mean_scale, P_scale, B_scale, L_scale mean_target, shape_scale, X_pro]...
    = shape_model(X_shape, mean(mass_areas), 0.98);
clear mean_X P_x B_x L_x mean_scale P_scale B_scale L_scale mean_target shape_scale;

figure('name', 'procrustes aligned large shapes'); hold on;
for ii = 1:51;
    plot(X_pro_l(ii,[1:500,1],1), X_pro_l(ii, [501:end,501]), colours(rem(ii,7)+1));
end    
%%
%
% Look at variances and mean errors from scale models

[mean(small_model.L_c) mean(large_model.L_c) mean(u_model.L_c); mean(r50_model.L_c) mean(r51_model.L_c) 0]

[mean(small_model.L_shape) mean(large_model.L_shape) mean(u_model.L_shape); mean(r50_model.L_shape) mean(r51_model.L_shape) 0]

[mean(small_model.L_tex) mean(large_model.L_tex) mean(u_model.L_tex); mean(r50_model.L_tex) mean(r51_model.L_tex) 0]

[mean(small_model.L_scale) mean(large_model.L_scale) mean(u_model.L_scale); mean(r50_model.L_scale) mean(r51_model.L_scale) 0]

[mean(er_small.ind.total) mean(er_large.ind.total) mean(er_u.ind.total); mean(er_r50.ind.total) mean(er_r51.ind.total) 0]

[mean(er_small.com.total) mean(er_large.com.total) mean(er_u.com.total); mean(er_r50.com.total) mean(er_r51.com.total) 0]

[mean(er_small.ind.shape) mean(er_large.ind.shape) mean(er_u.ind.shape); mean(er_r50.ind.shape) mean(er_r51.ind.shape) 0]

[mean(er_small.ind.tex) mean(er_large.ind.tex) mean(er_u.ind.tex); mean(er_r50.ind.tex) mean(er_r51.ind.tex) 0]

[mean(er_small.com.scale) mean(er_large.com.scale) mean(er_u.com.scale); mean(er_r50.com.scale) mean(er_r51.com.scale) 0]
%%
mo = mass_s.mass_outline;
m_diff = diff(mo);
mo = mo([true; m_diff(:,1) | m_diff(:,2)],:);
m_diff(~m_diff(:,1) & ~m_diff(:,2),:) = [];

m_lengths = [0; cumsum(sqrt(m_diff(:,1).^2 + m_diff(:,2).^2))];
%clear m_diff;

i_shape = interp1(m_lengths, mo, linspace(0, m_lengths(end),500));
%clear m_lengths;
    