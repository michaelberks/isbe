%%
for ii = 1:20
    go = 1; jj = 0;
    while go
        jj = jj+1;
        go = ~strcmp(m_files(jj).name, r_files(ii).name);
    end
    idx(ii) = jj;
end
pack
%%
MA = zeros(179,1);
for ii = 1:179
    load(['C:\isbe\dev\masses\', m_files(ii).name]);
    MA(ii) = mass.mass_area;
end
%%
[er_s] = model_errors2('mass_model', m_files, 'masses\', weights_s); save weights\er_s er_s;
[er_c] = model_errors2('mass_model', m_files, 'masses\', weights_c); save weights\er_c er_c;
[er_o3] = model_errors2('mass_model', m_files, 'masses\', weights_o3); save weights\er_o3 er_o3;
pack;

[er3_s] = model_errors3(m_files, weights_s); save weights\er3_s er3_s;
[er3_c] = model_errors3(m_files, weights_c); save weights\er3_c er3_c;
[er3_o3] = model_errors3(m_files, weights_o3); save weights\er3_o3 er3_o3;

pack;

[er_z] = model_errors2('mass_model', m_files, 'masses\', weights_z); save weights\er_z er_z;
[er_o] = model_errors2('mass_model', m_files, 'masses\', weights_o); save weights\er_o er_o;
[er_o2] = model_errors2('mass_model', m_files, 'masses\', weights_o2); save weights\er_o2 er_o2;

pack;


[er3_z] = model_errors3(m_files, weights_z); save weights\er3_z er3_z;
[er3_o] = model_errors3(m_files, weights_o); save weights\er3_o er3_o;
[er3_o2] = model_errors3(m_files, weights_o2); save weights\er3_o2 er3_o2;
%%
sd = std(RMS_scale);
mm = mean(RMS_scale);

figure('WindowStyle', 'docked'); hold on;

%plot([0 30], [mean(mm)+2*mean(sd) mean(mm)+2*mean(sd)], 'm');
for ii = 1:2
    
    plot([ii ii], [mm(ii)+2*sd(ii) mm(ii)-2*sd(ii)], 'bx', 'markersize', 5);
    plot([ii ii], [mm(ii)+2*sd(ii) mm(ii)-2*sd(ii)], 'b:', 'linewidth', 1.5);
    %plot(ii*ones(179,1), RMS_shape(:,ii), 'b.');
    plot(ii, mm(ii), 'rx', 'MarkerSize', 10);
end

%plot(0:30, [mean(mm) mm mean(mm)], 'r');
plot([0 3], [mean(mm) mean(mm)], 'g');
%%
mm = mean(ers);
figure('WindowStyle', 'docked'); hold on;
for ii = 1:3
    
    %plot([ii ii], [mm(ii)+2*sd(ii) mm(ii)-2*sd(ii)], 'bx', 'markersize', 5);
    %plot([ii ii], [mm(ii)+2*sd(ii) mm(ii)-2*sd(ii)], 'b:', 'linewidth', 1.5);
    plot(ii*ones(179,1), ers(:,ii), 'b.');
    plot(ii, mm(ii), 'rx', 'MarkerSize', 10);
end
%%
hold on;
boxplot(ers);
%plot(0:30, [mean(mm) mm mean(mm)], 'r');
%plot([0 30], [mean(mm) mean(mm)], 'g');
%%
figure('WindowStyle', 'docked'); hold on;
plot([1 1], [mean(er3_c.combined)+2*std(er3_c.combined)...
    mean(er3_c.combined)-2*std(er3_c.combined)], 'bx', 'markersize', 5);
plot([1 1], [mean(er3_c.combined)+2*std(er3_c.combined)...
    mean(er3_c.combined)-2*std(er3_c.combined)], 'b:', 'linewidth', 1.5);