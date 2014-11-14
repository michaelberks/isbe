generate_mass_AM(small_files, 'C:\isbe\dev\models\small_model', 500, 'C:\isbe\dev\masses\', 20000);
generate_mass_AM(large_files, 'C:\isbe\dev\models\large_model', 500, 'C:\isbe\dev\masses\', 20000);
generate_mass_AM(r_files50, 'C:\isbe\dev\models\r_model50', 500, 'C:\isbe\dev\masses\', 20000);
generate_mass_AM(r_files51, 'C:\isbe\dev\models\r_model51', 500, 'C:\isbe\dev\masses\', 20000);
%%
generate_mass_AM(r_files91, 'C:\isbe\dev\models\r_model91', 500, 'C:\isbe\dev\masses\', 20000);

%%
er_small = model_errors2('models\small_model', small_files, 'masses\'); save scale\er_small er_small;
%%
er_large = model_errors2('models\area_model', large_files, 'masses\'); save scale\er_large er_large;
%%
er_r50 = model_errors2('models\r_model50', r_files50, 'masses\'); save scale\er_r50 er_r50;
er_r51 = model_errors2('models\r_model51', r_files51, 'masses\'); save scale\er_r51 er_r51;
%%
generate_mass_AM(r_files138, 'C:\isbe\dev\models\r_model91', 500, 'C:\isbe\dev\masses\', 20000);
%%
generate_mass_AM(area_files, 'C:\isbe\dev\models\area_model', 500, 'C:\isbe\dev\masses\', 20000);
er_area = model_errors2('models\area_model', area_files, 'masses\'); save scale\er_area er_area;
%%
idx1 = zeros(101, 1);
for ii = 1:101
    k = randsample(length(unique_masses(ii).idx), 1);
    idx1(ii) = unique_masses(ii).idx(k);
end
u_files1 = m_files(idx1);
save C:\isbe\dev\files\u_files u_files1 idx1
%%
for ii = 1:20
    figure('WindowStyle', 'docked'); hold on;
    plot(X_shape(ii,1:500), X_shape(ii,501:1000), 'b');
    plot(ss(ii)*X_shape(ii,1:500), ss(ii)*X_shape(ii,501:1000), 'g');
    plot(mean_shape(1:500), mean_shape(501:1000), 'r');
end
%%
size_shape_vec = 500;
X_shape = zeros(101, 1000);
for ii = 1:101
    load(['C:\isbe\dev\masses\', u_files1(ii).name]);
    
    shape_vec = mass.mass_outline;
    idx = round(linspace(1, length(shape_vec(:,1)), size_shape_vec+1));
    idx = idx(1:size_shape_vec); %ensures first point is not equal to the last point!
    X_shape(ii,:) = [shape_vec(idx,1)', shape_vec(idx,2)'];
    clear mass;
end
%%
idx_small = zeros(91,1);
for ii = 1:91
    go = 1; jj = 0;
    while go
        jj = jj+1;
        go = ~strcmp(m_files(jj).name, small_files(ii).name);
    end
    idx_small(ii) = jj;
end
%%
mass_area = zeros(179,1);
for ii = 1:179
    load(['C:\isbe\dev\masses\', m_files(ii).name]);
    mass_area(ii) = mass.mass_area;
    clear mass;
end
%%
mm = mean(RMS_shape);
ss = std(RMS_shape);
figure('WindowStyle', 'docked'); hold on;
for ii = 1:29
    plot([ii ii], [mm(ii)+ss(ii), mm(ii)-ss(ii)], 'b:', 'LineWidth', 1.5);
    plot([ii ii], [mm(ii)+ss(ii), mm(ii)-ss(ii)], 'bx');
    plot([ii ii], [mm(ii), mm(ii)], 'rx', 'MarkerSize', 10);
end
plot([0, 30], [mean(mm), mean(mm)], 'g', 'LineWidth', 1.5);
%%
hold on;
plot([0.75, 1.25], [mean(ers(:,1)), mean(ers(:,1))], 'g', 'LineWidth', 1.5);
plot([1.75, 2.25], [mean(ers(:,2)), mean(ers(:,2))], 'g', 'LineWidth', 1.5);
plot([2.75, 3.25], [mean(ers(:,3)), mean(ers(:,3))], 'g', 'LineWidth', 1.5);
plot([3.75, 4.25], [mean(ers(:,4)), mean(ers(:,4))], 'g', 'LineWidth', 1.5);

plot([3.5, 4], [16 16], 'g', 'LineWidth', 1.5);
plot([3.5, 4], [15 15], 'r', 'LineWidth', 1.0);
%%
figure('WindowStyle', 'docked');
subplot(2,2,1);
bar(mean(abs(model_z.combined_data)), 'r');
subplot(2,2,2);
bar(mean(abs(model_s.combined_data)), 'g');
subplot(2,2,3);
bar(mean(abs(model_c.combined_data)), 'b');
subplot(2,2,4);
bar(mean(abs(model_o.combined_data)), 'c');
%%
figure('WindowStyle', 'docked');
subplot(2,2,1);
bar(mean(abs(l.mass_model.combined_data)), 'r');
subplot(2,2,2);
bar(mean(abs(s.mass_model.combined_data)), 'g');
subplot(2,2,3);
bar(mean(abs(m.mass_model.combined_data)), 'b');
subplot(2,2,4);
bar(mean(abs(mass_model.combined_data)), 'c');