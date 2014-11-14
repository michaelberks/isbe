% Work 12th March 2009

% Script transferring mass background structures into target regions

%Load in target region
args.TargetImage = imread('C:\isbe\dev\background\images\normal1024\o04_003LCC_1024_2541_743.bmp');

%Load in a mass background
load C:\isbe\dev\files\bg_files.mat
args.Mass = u_load(['C:\isbe\dev\masses\', bg_files(15).name]);

%Generate modified target region and display results
[modified_region modified_dt modified_c] = mb_transfer_mass_dual_tree(args);
figure('Name', 'Target region'); imagesc(args.TargetImage); axis image; colormap(gray(256));
figure('Name', 'Modified target region - no dilation'); imagesc(modified_region); axis image; colormap(gray(256));
hold on; plot(modified_c(1), modified_c(2), 'rx');
figure('Name', 'Difference map'); imagesc(modified_region - double(args.TargetImage)); axis image; colormap(gray(256));

%Now try using a larger mass region
args.MassDilate = 100;

[modified_region2 modified_dt2] = mb_transfer_mass_dual_tree(args);
figure('Name', 'Modified target region - dilation = 100'); imagesc(modified_region2); axis image; colormap(gray(256));
hold on; plot(modified_c(1), modified_c(2), 'rx');
figure('Name', 'Difference map'); imagesc(modified_region2 - double(args.TargetImage)); axis image; colormap(gray(256));
%%
%Now lets do this randomly for a load of images
load C:\isbe\dev\files\u_files.mat
target_list = dir('C:\isbe\dev\background\images\normal1024\*.bmp');
%%
args.MassDilate = 50;
bg_idx = randsample(1:length(target_list), 10);
for ii = 1:101
    args.Mass = u_load(['C:\isbe\dev\masses\', u_files1(ii).name]);
    if all(size(args.Mass.mass_ROI) <= [1024 1024])
        for jj = 1:10
            args.TargetImage = imread(['C:\isbe\dev\background\images\normal1024\', target_list(bg_idx(jj)).name]);
            args.Location = ([1024 1024] - fliplr(size(args.Mass.mass_ROI))) / 2;
            [modified_region modified_dt modified_c] = mb_transfer_mass_dual_tree(args);
            target_region = args.TargetImage;
            save(['C:\isbe\dev\background\modified_regions\mass', zerostr(ii,3), '_bg', zerostr(bg_idx(jj), 3)],...
                'modified_region', 'modified_dt', 'modified_c', 'target_region');
        end
    end
end

%%
region_list = dir('C:\isbe\dev\background\modified_regions\*.mat');
for ii = 401:410
    load(['C:\isbe\dev\background\modified_regions\', region_list(ii).name]);
    figure('Name', ['Modified target region - dilation = 50 ', region_list(ii).name]); 
    imagesc(modified_region); axis image; colormap(gray(256));
    hold on;
    plot(modified_c(1), modified_c(2), 'rx');
end

%%
% Now generate some new masses in the target zone
for ii = 61:101
    load(['C:\isbe\dev\background\modified_regions\', region_list(ii).name]);
    args.TargetRegion = modified_region;
    args.TargetCentre = modified_c;
    load C:\isbe\dev\mass_model\models\model_w500_50K.mat
    args.MassModel = mass_model;
    args.Plot = 1;
    [mass] = sample_new_mass_in_region(args);
    save(['C:\isbe\dev\background\new_masses\mass', zerostr(ii,3)], 'mass');
end
%%
% Now generate some masses but use the first
for ii = 1:length(region_list);
    load(['C:\isbe\dev\background\modified_regions\', region_list(ii).name]);
    mass_num = str2num(region_list(ii).name(end-6:end-4));
    args.TargetRegion = modified_region;
    args.TargetCentre = modified_c;
    args.ConditionMode = mass_num;
    load C:\isbe\dev\mass_model\models\model_w500_50K.mat
    args.MassModel = mass_model;
    args.Plot = 0;
    [mass] = sample_new_mass_in_region(args);
    save(['C:\isbe\dev\background\new_masses\mass', zerostr(ii,3)], 'mass');
end
%%
% These masses should have the same first combined mode as the real mass
% from which the background mass was sampled - yet they don't seem to
% matcht he background. Why? Is it just rotation (which doesn't match yet)
% or is there more? 
% Also: Does having a fixed dilation ring mean small masses are all
% increased to a generic circle? Should we dialte relative to mass size? Or
% better still mark new regions...
%
%%
% So lets build some new masses, directly in the real mass regions
for ii = 1:length(region_list);
    load(['C:\isbe\dev\background\modified_regions\', region_list(ii).name]);
    mass_num = str2num(region_list(ii).name(end-6:end-4));
    mass_real = u_load(['C:\isbe\dev\masses\', u_files1(mass_num).name]);
    args.TargetRegion = double(mass_real.mass_ROI) - mass_real.subtract_ROI;
    args.TargetCentre = mass_real.mass_centroid;
    args.ConditionMode = mass_num;
    load C:\isbe\dev\mass_model\models\model_w500_50K.mat
    args.MassModel = mass_model;
    args.Plot = 0;
    [mass] = sample_new_mass_in_region(args);
    save(['C:\isbe\dev\background\new_masses\real_copy5\mass', zerostr(ii,3)], 'mass');
end
%%
for ii = 1:30;
    mass_new = u_load(['C:\isbe\dev\background\new_masses\conditioned_scale1\mass', zerostr(ii,3)]);
    mass_num = str2num(region_list(ii).name(end-6:end-4));
    mass_real = u_load(['C:\isbe\dev\masses\', u_files1(mass_num).name]);
    cmax = max([mass_new.subtract_ROI(:); mass_real.subtract_ROI(:)]);
    figure; 
    subplot(1,3,1); imagesc(uint8(mass_new.subtract_ROI)); axis image; colormap(gray(256)); caxis([0 cmax]);
    subplot(1,3,2); imagesc(uint8(mass_real.subtract_ROI)); axis image; colormap(gray(256)); caxis([0 cmax]);
    subplot(1,3,3); imagesc(uint8(mass_new.mass_ROI)); axis image; colormap(gray(256));
end
%%
clear args
for ii = 1:30%length(region_list);
    load(['C:\isbe\dev\background\modified_regions\', region_list(ii).name]);
    mass_num = str2num(region_list(ii).name(end-6:end-4));
    mass_real = u_load(['C:\isbe\dev\masses\', u_files1(mass_num).name]);
    args.TargetRegion = double(mass_real.mass_ROI) - mass_real.subtract_ROI;
    args.TargetCentre = mass_real.mass_centroid;
    args.ScaleParameter = mass_model.B_scale(mass_num);
    args.ConditionMode = mass_num;
    args.ConditionLevel = 1;
    load C:\isbe\dev\mass_model\models\model_w500_50K.mat
    args.MassModel = generate_scale_model(mass_model);
    args.Plot = 0;
    [mass] = sample_scale_mass_in_region(args);
    save(['C:\isbe\dev\background\new_masses\conditioned_scale1\mass', zerostr(ii,3)], 'mass');
end
%%
%--------------------------------------------------------------------------
%For each new mass
for ii = 1:30;
    mass_new = u_load(['C:\isbe\dev\background\new_masses\mass', zerostr(ii,3)]);
    mass_num = str2num(region_list(ii).name(end-6:end-4));
    mass_real = u_load(['C:\isbe\dev\masses\', u_files1(mass_num).name]);
    cmax = max([mass_new.subtract_ROI(:); mass_real.subtract_ROI(:)]);
    figure; 
    subplot(1,3,1); imagesc(uint8(mass_new.subtract_ROI)); axis image; colormap(gray(256)); caxis([0 cmax]);
    subplot(1,3,2); imagesc(uint8(mass_real.subtract_ROI)); axis image; colormap(gray(256)); caxis([0 cmax]);
    subplot(1,3,3); imagesc(uint8(mass_new.mass_ROI)); axis image; colormap(gray(256));
end

%%
% Tomorrow - condition model based on scale

%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%%
% Lets use all the masses from the dataset - is independence now really
% that important? Yes! It is, for the same reason building the model from
% independent masses is a good idea
% However the masses we've not been using in our models have outline
% centered on [0 0] instead of the centroid - so change this
old_mass = dir('C:\isbe\dev\masses\temp\*.mat');
for ii = 11:length(old_mass);
    mass = u_load(['C:\isbe\dev\masses\temp\', old_mass(ii).name]);
    mass.mass_outline = mass.mass_outline + repmat(mass.mass_centroid, size(mass.mass_outline, 1), 1);
    %figure; imagesc(mass.mass_ROI); axis image; colormap(gray(256));
    %hold on;
    %plot(mass.mass_outline(:,1), mass.mass_outline(:,2), 'c');
    save(['C:\isbe\dev\masses\temp\', old_mass(ii).name], 'mass');
end
%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%Can also try again with half-size regions - though I think we should give up on this now
args.TargetImage = imread('C:\isbe\dev\background\images\normal512\o04_003LCC_1024_2541_743.bmp');
args.CoarsestLevel = 4;
%Need to change some code in function (uncomment the 0.5 imresize steps)
[modified_region3 modified_dt3] = mb_transfer_mass_dual_tree(args);