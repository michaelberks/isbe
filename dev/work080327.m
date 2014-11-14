%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thursday March 27th
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Generate new masses
clear
new_mass_args.MassModel = u_load('C:\isbe\dev\models\model_w500_50K.mat');
new_mass_args.NumberOfMasses = 101;
new_mass_args.SavePath = 'C:\isbe\dev\new_masses\masses080327\';
sample_new_masses(new_mass_args);

%%
%Now compute summary statistics for the sets
real_files = u_load('C:\isbe\dev\files\u_files.mat');
syn_files = dir('C:\isbe\dev\new_masses\masses080327\new_mass*');

[S_real] = mass_statistics('MassNames', real_files, 'MassPath', 'C:\isbe\dev\masses\');
[S_syn] = mass_statistics('MassNames', syn_files, 'MassPath', 'C:\isbe\dev\new_masses\masses080327\');
save C:\isbe\dev\new_masses\mass_stats S_*
%%
%Display the statistics and do some tests
display('Means of summary statistics:');
display(['Area: Real = ', num2str(mean(S_real.AreaMM)), ' Syn = ', num2str(mean(S_syn.AreaMM))]);
display(['Perimeter: Real = ', num2str(mean(S_real.PerimeterMM)), ' Syn = ', num2str(mean(S_syn.PerimeterMM))]);
display(['Texture mean: Real = ', num2str(mean(S_real.TextureMean)), ' Syn = ', num2str(mean(S_syn.TextureMean))]);
display(['Texture max: Real = ', num2str(mean(S_real.TextureMax)), ' Syn = ', num2str(mean(S_syn.TextureMax))]);
display(['Texture variance: Real = ', num2str(mean(S_real.TextureVar)), ' Syn = ', num2str(mean(S_syn.TextureVar))]);
%Display the statistics and do some tests
display('Standard errors of summary statistics:');
display(['Area: Real = ', num2str(std(S_real.AreaMM)/sqrt(101)), ' Syn = ', num2str(std(S_syn.AreaMM)/sqrt(101))]);
display(['Perimeter: Real = ', num2str(std(S_real.PerimeterMM)/sqrt(101)), ' Syn = ', num2str(std(S_syn.PerimeterMM)/sqrt(101))]);
display(['Texture mean: Real = ', num2str(std(S_real.TextureMean)/sqrt(101)), ' Syn = ', num2str(std(S_syn.TextureMean)/sqrt(101))]);
display(['Texture max: Real = ', num2str(std(S_real.TextureMax)/sqrt(101)), ' Syn = ', num2str(std(S_syn.TextureMax)/sqrt(101))]);
display(['Texture variance: Real = ', num2str(std(S_real.TextureVar)/sqrt(101)), ' Syn = ', num2str(std(S_syn.TextureVar)/sqrt(101))]);
%%
display('Hypothesis tests (unpaired t-test) for stats');
[h,p] = ttest2(S_real.Area, S_syn.Area);
if h
    display(['Area: Significant difference, p = ', num2str(p)]);
else
    display(['Area: No significant difference, p = ', num2str(p)]);
end

[h,p] = ttest2(S_real.Perimeter, S_syn.Perimeter);
if h
    display(['Perimeter: Significant difference, p = ', num2str(p)]);
else
    display(['Perimeter: No significant difference, p = ', num2str(p)]);
end

[h,p] = ttest2(S_real.TextureMean, S_syn.TextureMean);
if h
    display(['Texture mean: Significant difference, p = ', num2str(p)]);
else
    display(['Texture mean: No significant difference, p = ', num2str(p)]);
end

[h,p] = ttest2(S_real.TextureMax, S_syn.TextureMax);
if h
    display(['Texture max: Significant difference, p = ', num2str(p)]);
else
    display(['Texture max: No significant difference, p = ', num2str(p)]);
end

[h,p] = ttest2(S_real.TextureVar, S_syn.TextureVar);
if h
    display(['Texture var: Significant difference, p = ', num2str(p)]);
else
    display(['Texture var: No significant difference, p = ', num2str(p)]);
end
%%
% Ok, we're happy the masses have statistics similar enough to real masses
% let's look at them then...
for ii = 61:80
    mass = u_load(['C:\isbe\dev\new_masses\masses080327\', syn_files(ii).name]);
    figure('Name', ['Mass ', num2str(ii)]); 
    imagesc(uint8(mass.subtract_ROI)); axis image; colormap(gray(256)); colorbar;
    clear mass;
end
% Masses that look (really) great:
% [71, lots more than this; just can't choose yet!]
% Masses that look crap!
% [32, 46, 49, 53, 78, 83, 88]

%%
%find each mass's buddy
real_buddy = zeros(101,1);
for ii = 1:101
    [dummy idx] = min(abs(S_real.Area - S_syn.Area(ii)));
    real_buddy(ii) = idx;
    clear idx dummy
end

%%
for ii = 1:20
    syn_mass = u_load(['C:\isbe\dev\new_masses\masses080327\', syn_files(ii).name]);
    real_mass = u_load(['C:\isbe\dev\masses\', real_files(real_buddy(ii)).name]);
    
    cmax = max([S_real.TextureMax(real_buddy(ii)), S_syn.TextureMax(ii)]);
    
    figure('Name', ['Mass ', num2str(ii)]); 
    subplot(1,2,1);
    imagesc(uint8(real_mass.subtract_ROI)); axis image; colormap(gray(256)); caxis([0 cmax]);
    subplot(1,2,2);
    imagesc(uint8(syn_mass.subtract_ROI)); axis image; colormap(gray(256)); caxis([0 cmax]);
    clear real_mass syn_mass;
end
%%
%Look at the 10 synthetic masses with highest mean intensity
for ii = 1:10
    idx = b_idx(ii);
    syn_mass = u_load(['C:\isbe\dev\new_masses\masses080327\', syn_files(idx).name]);
    figure('Name', ['Mass ', num2str(idx)]);
    imagesc(uint8(syn_mass.subtract_ROI)); axis image; colormap(gray(256)); caxis([0 80]);
end
%%
% Look at the 53 nice background - try all at once on this comp
for ii = 1:53
    real_mass = u_load(['C:\isbe\dev\masses\', bg_files(ii).name]);
    figure('Name', ['Mass ', bg_files(ii).name]);
    imagesc(double(real_mass.mass_ROI) - real_mass.subtract_ROI); axis image; colormap(gray(256));
end
%%
for ii = 6:15
    load(['C:\isbe\dev\background\syn\mass_2_2levels\conditioned_synthesis_circle_6_', zerostr(ii,3), '.mat']);
    figure; imagesc(synthesised_image); colormap(gray(256)); axis image;
    clear synthesised_image pyramid cluster_image;
end
%%
% Failed synthesis: 8, 10, 12, 13, 14 - which sub-band(s) failed?

% Synthesis 8:
load('C:\isbe\dev\background\syn\mass_2_2levels\conditioned_synthesis_circle_6_008.mat');
for lev = 2:4
    for ori = 1:5
        figure; imagesc(pyramid{lev, ori}); colormap(jet(256)); axis image;
    end
end
%failure in {2,3}, {2,4}, {3,3}, {3,4}
%%
% Synthesis 10:
load('C:\isbe\dev\background\syn\mass_2_2levels\conditioned_synthesis_circle_6_010.mat');
for lev = 2:6
    for ori = 1:5
        figure; imagesc(pyramid{lev, ori}); colormap(jet(256)); axis image;
    end
end
%failure in {2,3}, {2,4}
%%
% Synthesis 12:
load('C:\isbe\dev\background\syn\mass_2_2levels\conditioned_synthesis_circle_6_012.mat');
for lev = 2:4
    for ori = 1:5
        figure; imagesc(pyramid{lev, ori}); colormap(jet(256)); axis image;
    end
end
%failure in {2,4}
%%
% Synthesis 13:
load('C:\isbe\dev\background\syn\mass_2_2levels\conditioned_synthesis_circle_6_013.mat');
for lev = 2:6
    for ori = 1:5
        figure; imagesc(pyramid{lev, ori}); colormap(jet(256)); axis image;
    end
end
%failure in {2,3}, {2,4}
%%
% Synthesis 14:
load('C:\isbe\dev\background\syn\mass_2_2levels\conditioned_synthesis_circle_6_014.mat');
for lev = 2:4
    for ori = 1:5
        figure; imagesc(pyramid{lev, ori}); colormap(jet(256)); axis image;
    end
end
%failure in {2,3}, {2,4}
%%
for ii = 1:15
    load(['C:\isbe\dev\background\syn\mass_2_2levels\conditioned_synthesis', zerostr(ii,3), '.mat']);
    i1 = double(imread(['C:\isbe\dev\background\images\normal_2\normal', zerostr(ii,3), '.bmp']));
    figure; 
    subplot(1,2,1); imagesc(synthesised_image); colormap(gray(256)); axis image;
    subplot(1,2,2); imagesc(i1 - synthesised_image); colormap(gray(256)); axis image;
    clear synthesised_image pyramid cluster_image i1;
end
%%
%we've added masses to subtracted and normal regions. Now compare them and
%pick suitable example for paper

mass6_regions = dir('C:\isbe\dev\new_masses\regions\mass06*');
mass44_regions = dir('C:\isbe\dev\new_masses\regions\mass44*');
mass64_regions = dir('C:\isbe\dev\new_masses\regions\mass64*');
%%
%mass 6
for ii = 1:length(mass6_regions)
    load(['C:\isbe\dev\new_masses\regions\', mass6_regions(ii).name]);
    figure('Name', ['Region ', mass6_regions(ii).name]) ; 
    imagesc(mass.background_ROI + mass.subtract_ROI); axis image; colormap(gray(256));
    clear mass;
end
%%
%mass 44
for ii = 1:length(mass44_regions)
    load(['C:\isbe\dev\new_masses\regions\', mass44_regions(ii).name]);
    figure('Name', ['Region ', mass44_regions(ii).name]) ; 
    imagesc(mass.background_ROI + mass.subtract_ROI); axis image; colormap(gray(256));
    clear mass;
end
%%
%mass 64
for ii = 1:length(mass64_regions)
    load(['C:\isbe\dev\new_masses\regions\', mass64_regions(ii).name]);
    figure('Name', ['Region ', mass64_regions(ii).name]) ; 
    imagesc(mass.background_ROI + mass.subtract_ROI); axis image; colormap(gray(256));
    clear mass;
end
%%
% I like mass 64, bg 06, norm 15
mass_bg = u_load('C:\isbe\dev\new_masses\regions\mass64_bg06.mat');
mass_norm = u_load('C:\isbe\dev\new_masses\regions\mass64_norm15.mat');
figure; imagesc(mass_bg.background_ROI + mass_bg.subtract_ROI); axis image; colormap(gray(256)); %caxis([80 225]);
figure; imagesc(mass_bg.background_ROI); axis image; colormap(gray(256)); %caxis([80 225]);
figure; imagesc(mass_norm.background_ROI + mass_norm.subtract_ROI); axis image; colormap(gray(256)); %caxis([80 225]);
figure; imagesc(mass_norm.background_ROI); axis image; colormap(gray(256)); %caxis([80 225]);