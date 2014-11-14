%%
for ii = 1:length(files)
    mass = load(['C:\isbe\mammograms\new_CAD\annotations\', files(ii).name]);
    name = ['C:\isbe\mammograms\new_CAD\BMP_2004\o', ...
            files(ii).name(3:11), '.bmp'];
    mass.name = name;
    save(['C:\isbe\mammograms\new_CAD\annotations\', files(ii).name],...
            'mass');
    clear mass name;
end
clear ii
%%
for ii = 1:15
    load(['C:\isbe\mammograms\new_CAD\annotations\', sub_masses(ii).name]);
    f(ii) = figure('WindowStyle', 'docked', 'name', sub_masses(ii).name);
    imagesc(mass.mass_ROI); colormap(jet(256)); hold on; axis image;
    p(ii) = plot(mass.mass_outline([1:end,1],1),...
        mass.mass_outline([1:end,1],2), 'r');
end

%%
for ii = 179:265
    i1 = imread(files(ii).name);
    f_name = [files(ii).name(1:11), 'jpg'];
    imwrite(i1, f_name);
    clear i1 f_name;
end
%%
sub_masses(1).name = 'an04_009RML.mat';
sub_masses(2).name = 'an04_018LCC.mat';
sub_masses(3).name = 'an04_031LML.mat';
sub_masses(4).name = 'an04_036RML.mat';
sub_masses(5).name = 'an04_038LCC.mat';
sub_masses(6).name = 'an04_039RCC.mat';
sub_masses(7).name = 'an04_040RMLb.mat';
sub_masses(8).name = 'an04_045LML.mat';
sub_masses(9).name = 'an04_046RCC.mat';
sub_masses(10).name = 'an04_057LML.mat';
sub_masses(11).name = 'an04_059RCC.mat';
sub_masses(12).name = 'an04_074RCC.mat';
sub_masses(13).name = 'an04_075RCC.mat';
sub_masses(14).name = 'an04_081LMLa.mat';
sub_masses(15).name = 'an04_131RML.mat';
save C:\isbe\dev\subtraction\sub_masses sub_masses;

%%
cases(1).name = files(1).name(6:9);
for ii = 2:185
    if ~strncmp(files(ii-1).name, files(ii).name, 9)
        cases(end+1).name = files(ii).name(6:9);
    end
end
for ii = 1:length(cases)
    display(cases(ii).name);
end
    
%%
bg = imread('C:\isbe\dev\subtraction\bg\dense.bmp');
for ii = 1:15
    if ~rem(ii-1, 12)
        figure('WindowStyle', 'docked');
    end
    load(['C:\isbe\mammograms\new_CAD\annotations\', sub_masses(ii).name]);
    subplot(3,4, rem(ii-1, 12)+1);
    cx = min(mass.mass_outline(:,1)) + 0.5*(max(mass.mass_outline(:,1))...
        - min(mass.mass_outline(:,1)));
    cy = min(mass.mass_outline(:,2)) + 0.5*(max(mass.mass_outline(:,2))...
        - min(mass.mass_outline(:,2)));
    image(bg); axis image; colormap(gray(256)); hold on;
    plot(mass.mass_outline([1:end,1],1) + 750 - cx,...
        mass.mass_outline([1:end,1],2) + 750 - cy, 'r');
    title(sub_masses(ii).name);
    clear mass;
end
