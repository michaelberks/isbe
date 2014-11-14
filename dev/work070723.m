%%
% Left

load C:\isbe\dev\annotations\an04_137LCC
d = 20;
mass.mass_ROI(:,1:d) = [];
mass.C1 = mass.C1 + d;
mass.mass_outline(:,1) = mass.mass_outline(:,1) - d;
for ii = 1:length(mass.mass_spicules)
    mass.mass_spicules(ii).outline(:,1) = ...
        mass.mass_spicules(ii).outline(:,1) - d;
end

save C:\isbe\dev\temp\an04_137LCC

%%
% Right

load C:\isbe\dev\annotations\an04_051RCCb
d = 360;
mass.mass_ROI(:,d+1:end) = [];
mass.C2 = mass.C1 + d - 1;

save C:\isbe\dev\temp\an04_051RCCb
clear mass d;
%%
% Plot
load C:\isbe\dev\annotations\an04_141RCCb

figure('WindowStyle', 'docked');
subplot(1,2,1)
imagesc(mass.mass_ROI); axis image; colormap(gray(256)); hold on;
plot(mass.mass_outline(:,1), mass.mass_outline(:,2), 'r:');
for ii = 1:length(mass.mass_spicules)
    plot(mass.mass_spicules(ii).outline(:,1),...
        mass.mass_spicules(ii).outline(:,2), 'r:');
end

load C:\isbe\dev\temp\an04_141RCCb

subplot(1,2,2)
imagesc(mass.mass_ROI); axis image; colormap(gray(256)); hold on;
plot(mass.mass_outline(:,1), mass.mass_outline(:,2), 'r:');
for ii = 1:length(mass.mass_spicules)
    plot(mass.mass_spicules(ii).outline(:,1),...
        mass.mass_spicules(ii).outline(:,2), 'r:');
end
clear mass;
%%
d = 360;
load C:\isbe\dev\temp\an04_051RCCb
figure;
imagesc(mass.mass_ROI); axis image; colormap(gray(256)); hold on;
plot([d d], [1 size(mass.mass_ROI, 1)], 'r:');
plot(mass.mass_outline(:,1), mass.mass_outline(:,2), 'gx');
plot(mass.mass_outline(1,1), mass.mass_outline(1,2), 'ro');
plot(mass.mass_outline(2,1), mass.mass_outline(2,2), 'bo');
%%
d = 360;
load C:\isbe\dev\temp\an04_051RCCb

idx = mass.mass_outline(:,1) > d - 1;
mass.mass_outline(idx, 1) = d - 1;

figure;
imagesc(mass.mass_ROI); axis image; colormap(gray(256)); hold on;
plot(mass.mass_outline(:,1), mass.mass_outline(:,2), 'r:');

save C:\isbe\dev\temp\an04_051RCCb
%%
edge = [12 13 16 25 26 27 28 31 61 67 69 86 94 98 99 106 117 118 119 120 123 125 128 165 178 181];
e_files = files(edge);
for ii = 1:length(e_files);
    display([e_files(ii).name]);
end
%%
for ii = [4 6 8 10 15 16 21 24 26]
    load(['C:\isbe\dev\temp\', e_files(ii).name]);
    figure('WindowStyle', 'docked');
    subplot(1,2,1)
    imagesc(mass.mass_spline); axis image; colormap(jet(256)); hold on;
    
    subplot(1,2,2)
    load(['C:\isbe\dev\annotations\', e_files(ii).name]);
    imagesc(mass.mass_spline); axis image; colormap(jet(256)); hold on;
end
%%
subtract_mass(files, 'C:\isbe\dev\annotations\', 50, 10, 10, 18, 0, 'biharmTPS', 1);
subtract_mass_it(files, 'C:\isbe\dev\annotations\', 50, 10, 10, 18, 0, 'biharmTPS', 1);
%%

for ii = 1:179
    load(['C:\isbe\dev\annotations\', files(ii).name]);
    [mv vv] = mass_halo(mass, 24, 50, 20, max(0, 11-ii));
    
    load(['C:\isbe\dev\masses\', m_files(ii).name]);
    mass.mass_halo = [mv; vv];
    save(['C:\isbe\dev\masses\', m_files(ii).name], 'mass');
    clear mass;
end
%%
for ii = 1:179
    load(['C:\isbe\dev\masses\', m_files(ii).name]);
    if sum(isnan(mass.mass_halo)); break; end
    if sum(isinf(mass.mass_halo)); break; end
    clear mass;
end