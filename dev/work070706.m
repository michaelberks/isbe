for ii = 81:102
    temp = load(['C:\isbe\mammograms\new_CAD\rumana\mb\', a_files(ii).name]);
    mass_r = temp.mass; clear temp;
    
    temp = load(['C:\isbe\mammograms\new_CAD\rumana\annotations\', a_files(ii).name]);
    mass_a = temp.mass; clear temp;
    
    temp = load(['C:\isbe\mammograms\new_CAD\rumana\changed\', c_files(ii).name]);
    mass_c = temp.mass; clear temp;
    
    if (mass_c.R1 ~= mass_r.R1)
        mass_c.mass_outline(:,2) = mass_c.mass_outline(:,2)...
            + mass_c.R1 - mass_r.R1;
        %display(['y shifted in mass ', num2str(ii)])
    end
    if (mass_c.C1 ~= mass_r.C1)
        mass_c.mass_outline(:,1) = mass_c.mass_outline(:,1)...
            + mass_c.C1 - mass_r.C1;
        %display(['x shifted in mass ', num2str(ii)])
    end
    if (mass_a.R1 ~= mass_r.R1)
        mass_a.mass_outline(:,2) = mass_a.mass_outline(:,2)...
            + mass_a.R1 - mass_r.R1;
        %display(['y shifted in mass ', num2str(ii)])
    end
    if (mass_a.C1 ~= mass_r.C1)
        mass_a.mass_outline(:,1) = mass_a.mass_outline(:,1)...
            + mass_a.C1 - mass_r.C1;
        %display(['x shifted in mass ', num2str(ii)])
    end
    figure('WindowStyle', 'docked');
%     subplot(1,2,1)
    imagesc(mass_r.mass_ROI); colormap(gray(256)); axis('image'); hold('on');
    plot(mass_r.mass_outline(:,1), mass_r.mass_outline(:,2), 'r');
    
%     subplot(1,2,2)
%     imagesc(mass_c.mass_ROI); colormap(gray(256)); axis('image'); hold('on');
    plot(mass_c.mass_outline(:,1), mass_c.mass_outline(:,2), 'b:');
    plot(mass_a.mass_outline(:,1), mass_a.mass_outline(:,2), 'y:');
    title(c_files(ii).name);
    clear mass_c mass_r mass_a;
end
%%
for ii = 1:20
    load(['changed\', c_files(ii).name]);
    
    mo = mass.mass_outline'; clear mass;
    mo(:,end+1) = mo(:,1);
    
    t = 1:size(mo,2);
    tt = 1:0.2:size(mo,2);
    mos = spline(t, mo, tt);
    
    figure('WindowStyle', 'docked');
    plot(mos(1,:), mos(2,:));
    hold on;
    plot(mo(1,:), mo(2,:), 'r:');
end
%%
    temp = load(['C:\isbe\mammograms\new_CAD\rumana\annotations\an04_008RCCRR.mat']);
    mass_r = temp.mass; clear temp;
    
    temp = load('C:\isbe\mammograms\new_CAD\rumana\changed\an04_008RCC.mat');
    mass_c = temp.mass; clear temp;
    
    if (mass_c.R1 ~= mass_r.R1)
        mass_c.mass_outline(:,2) = mass_c.mass_outline(:,2)...
            + mass_c.R1 - mass_r.R1;
    end
    if (mass_c.C1 ~= mass_r.C1)
        mass_c.mass_outline(:,1) = mass_c.mass_outline(:,1)...
            + mass_c.C1 - mass_r.C1;
    end
    figure('WindowStyle', 'docked');
    imagesc(mass_r.mass_ROI); colormap(gray(256)); axis('image'); hold('on');
    plot(mass_r.mass_outline(:,1), mass_r.mass_outline(:,2), 'r');
    plot(mass_c.mass_outline(:,1), mass_c.mass_outline(:,2), 'b:');
    
    clear mass_c mass_r;
%%
i1 = imread('C:\isbe\mammograms\new_CAD\BMP_2004\o04_008RML.bmp');
figure('WindowStyle', 'docked');
imagesc(i1); axis('image'); colormap(gray(256)); hold('on');
clear i1;
load C:\isbe\mammograms\new_CAD\annotations\an04_008RML.mat;
plot(mass.mass_outline(:,1)+mass.C1, mass.mass_outline(:,2)+mass.R1, 'b:');
clear mass;
load C:\isbe\mammograms\new_CAD\annotations\an04_008RMLsmall.mat;
plot(mass.mass_outline(:,1)+mass.C1, mass.mass_outline(:,2)+mass.R1, 'r:');
clear mass;
%%
load C:\isbe\mammograms\new_CAD\rumana\mb\an04_040RCCaRR.mat;
plot(mass.mass_outline(:,1)+mass.C1, mass.mass_outline(:,2)+mass.R1, 'r:');
clear mass;