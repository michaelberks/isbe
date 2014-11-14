%%
for ii = 1:51
    file_out = ['C:\isbe\dev\models\loo_r51\model', zerostr(ii, 3)];
    one_out_files = r_files51;
    one_out_files(ii) = [];

    mass_model = generate_mass_AM(...
        one_out_files, file_out, 500, 'C:\isbe\dev\masses\', 20000);
end
[er_r51.com er_r51.ind] = model_errors_loo(r_files51, [], 'C:\isbe\dev\models\loo_r51\model');
save C:\isbe\dev\scale\er_r51 er_r51;
%%

load C:\isbe\dev\files\m_files.mat
%%
figure('WindowStyle', 'docked');
for ii = 1:20
    load(['C:\isbe\dev\masses\', m_files(ii).name]);
    subplot(4,5,ii)
    plot(mass.mass_outline([1:end,1],1), mass.mass_outline([1:end,1],2), 'k');
    title(m_files(ii).name);
    saveas(gcf, 'C:\isbe\dev\figures\masses1to20.eps');
end

figure('WindowStyle', 'docked');
for ii = 21:40
    load(['C:\isbe\dev\masses\', m_files(ii).name]);
    subplot(4,5,ii-20)
    plot(mass.mass_outline([1:end,1],1), mass.mass_outline([1:end,1],2), 'k');
    title(m_files(ii).name);
    saveas(gcf, 'C:\isbe\dev\figures\masses21to40.eps');
end

figure('WindowStyle', 'docked');
for ii = 41:60
    load(['C:\isbe\dev\masses\', m_files(ii).name]);
    subplot(4,5,ii-40)
    plot(mass.mass_outline([1:end,1],1), mass.mass_outline([1:end,1],2), 'k');
    title(m_files(ii).name);
    saveas(gcf, 'C:\isbe\dev\figures\masses41to60.eps');
end

figure('WindowStyle', 'docked');
for ii = 61:80
    load(['C:\isbe\dev\masses\', m_files(ii).name]);
    subplot(4,5,ii-60)
    plot(mass.mass_outline([1:end,1],1), mass.mass_outline([1:end,1],2), 'k');
    title(m_files(ii).name);
    saveas(gcf, 'C:\isbe\dev\figures\masses61to80.eps');
end

figure('WindowStyle', 'docked');
for ii = 81:100
    load(['C:\isbe\dev\masses\', m_files(ii).name]);
    subplot(4,5,ii-80)
    plot(mass.mass_outline([1:end,1],1), mass.mass_outline([1:end,1],2), 'k');
    title(m_files(ii).name);
    saveas(gcf, 'C:\isbe\dev\figures\masses81to100.eps');
end

figure('WindowStyle', 'docked');
for ii = 101:120
    load(['C:\isbe\dev\masses\', m_files(ii).name]);
    subplot(4,5,ii-100)
    plot(mass.mass_outline([1:end,1],1), mass.mass_outline([1:end,1],2), 'k');
    title(m_files(ii).name);
    saveas(gcf, 'C:\isbe\dev\figures\masses101to120.eps');
end

figure('WindowStyle', 'docked');
for ii = 121:140
    load(['C:\isbe\dev\masses\', m_files(ii).name]);
    subplot(4,5,ii-120)
    plot(mass.mass_outline([1:end,1],1), mass.mass_outline([1:end,1],2), 'k');
    title(m_files(ii).name);
    saveas(gcf, 'C:\isbe\dev\figures\masses121to140.eps');
end

figure('WindowStyle', 'docked');
for ii = 141:160
    load(['C:\isbe\dev\masses\', m_files(ii).name]);
    subplot(4,5,ii-140)
    plot(mass.mass_outline([1:end,1],1), mass.mass_outline([1:end,1],2), 'k');
    title(m_files(ii).name);
    saveas(gcf, 'C:\isbe\dev\figures\masses141to160.eps');
end

figure('WindowStyle', 'docked');
for ii = 161:179
    load(['C:\isbe\dev\masses\', m_files(ii).name]);
    subplot(4,5,ii-160)
    plot(mass.mass_outline([1:end,1],1), mass.mass_outline([1:end,1],2), 'k');
    title(m_files(ii).name);
    saveas(gcf, 'C:\isbe\dev\figures\masses161to179.eps');
end