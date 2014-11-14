% iwdm_poster _script
% Prepares the images to be shown on the poster for iwdm

%1) Save some colour images of a pyramid decomposition
load('C:\isbe\dev\background\pyramid5\normal512\o04_001LCC_1024_3427_865_pyramid.mat')
mkdir C:\isbe\dev\background\misc\pyr_ims\o04_001LCC
for level = 2:6
    for ori = 1:5
        write_im_from_colormap(pyr{level, ori}, ['C:\isbe\dev\background\misc\pyr_ims\o04_001LCC\subband_', num2str(level), '_', num2str(ori), '.bmp'], jet(256), [-5 5]);
    end
end
i1 = imread('C:\isbe\dev\background\images\normal512\o04_001LCC_1024_3427_865.bmp');
[gp p_sizes] = buildGpyr(double(i1), 5);
gp = mb_change_pyramid_form(gp, p_sizes, 'g');

for level = 1:5
    write_im_from_colormap(gp{level}, ['C:\isbe\dev\background\misc\pyr_ims\o04_001LCC\gp_', num2str(level), '.bmp'], gray(256), [0 256]);
end

%%
%2) Make synthetic masses and real masses for poster border

%make a directory for them
mkdir K:\isbe\conferences_and_symposia\iwdm2008\poster\figures
%%
%i) Select 6 subtracted mass backgrounds of sufficient size (756x756)

load C:\isbe\dev\files\bg_files.mat
load C:\isbe\dev\files\m_files.mat
load C:\isbe\dev\files\u_files.mat
%%
bg_to_add_mass_idx = [];
bg_regions = {};
for ii = 1:53
    mass = u_load(['C:\isbe\dev\masses\', bg_files(ii).name]);
    anno = u_load(['C:\isbe\dev\annotations\an', bg_files(ii).name(2:end)]);
    mammo = double(imread(anno.name));
    start_dims = round(mass.mass_centroid - 377);
    if all(start_dims>=1) && all(start_dims + 755 < min(size(mass.mass_ROI)))
        bg_to_add_mass_idx(end+1) = ii; %#ok
        bg_regions{end+1} = double(mass.mass_ROI(start_dims(2):start_dims(2)+755, start_dims(1):start_dims(1)+755))...
            - mass.subtract_ROI(start_dims(2):start_dims(2)+755, start_dims(1):start_dims(1)+755); %#ok
    elseif all(start_dims + min(anno.R1, anno.C1) >=1) && all(start_dims + max(anno.R1, anno.C1) + 755 < min(size(mammo)))
        bg_to_add_mass_idx(end+1) = ii; %#ok        
        if strcmpi(bg_files(ii).name(8), 'R')
            mammo(anno.R1:anno.R2, anno.C1:anno.C2) = fliplr(double(mass.mass_ROI) - mass.subtract_ROI);
        else
            mammo(anno.R1:anno.R2, anno.C1:anno.C2) = double(mass.mass_ROI) - mass.subtract_ROI;
        end
        bg_regions{end+1} = mammo(start_dims(2)+anno.R1-1:start_dims(2)+anno.R1+754, start_dims(1)+anno.C1-1:start_dims(1)+anno.C1+754); %#ok
    else
        display(['Mass ', num2str(ii), ' won''t fit']);        
    end
    clear mass anno mammo
end
%%
save K:\isbe\conferences_and_symposia\iwdm2008\poster\figures\bg_regions bg_regions
bg_mass_idx = randsample(bg_to_add_mass_idx, 6);
bg_mass_idx = sort(bg_mass_idx);

%%
%ii) Select 12 real masses of sufficient size excluding those used for
%synthetic mass backgrounds
idx_u2 = setdiff(idx_u1, idx_bg(bg_mass_idx));
u_files2 = m_files(idx_u2);

real_mass_idx = [];
mass_regions = {};
for ii = 1:length(u_files2)
    load(['C:\isbe\dev\masses\', u_files2(ii).name]);
    anno = u_load(['C:\isbe\dev\annotations\an', u_files2(ii).name(2:end)]);
    mammo = double(imread(anno.name));
    start_dims = round(mass.mass_centroid - 377);
    if all(start_dims>=1) && all(start_dims + 755 < min(size(mass.mass_ROI)))
       real_mass_idx(end+1) = ii; %#ok
       mass_regions{end+1} = double(mass.mass_ROI(start_dims(2):start_dims(2)+755, start_dims(1):start_dims(1)+755)); %#ok
    elseif all(start_dims + min(anno.R1, anno.C1) >=1) && all(start_dims + max(anno.R1, anno.C1) + 755 < min(size(mammo)))
        real_mass_idx(end+1) = ii; %#ok        
        if strcmpi(u_files2(ii).name(8), 'R')
            mammo(anno.R1:anno.R2, anno.C1:anno.C2) = fliplr(double(mass.mass_ROI));
        else
            mammo(anno.R1:anno.R2, anno.C1:anno.C2) = double(mass.mass_ROI);
        end
        mass_regions{end+1} = mammo(start_dims(2)+anno.R1-1:start_dims(2)+anno.R1+754, start_dims(1)+anno.C1-1:start_dims(1)+anno.C1+754); %#ok
    else
        display(['Mass ', num2str(ii), ' won''t fit']);        
    end
    clear mass anno mammo
end
save K:\isbe\conferences_and_symposia\iwdm2008\poster\figures\mass_regions mass_regions
real_mass_idx = sort(randsample(1:73, 13));
%%
for ii = 1:13    
    figure; imagesc(mass_regions{real_mass_idx(ii)}); axis image; caxis([10, 255]); colormap(gray(256));
    write_im_from_colormap(mass_regions{real_mass_idx(ii)}, ['K:\isbe\conferences_and_symposia\iwdm2008\poster\figures\real_mass', zerostr(ii,2), '.bmp'], gray(256), [10 255]);
end
%%
%iii) Select 6 normal regions and take the top 756x756 portion
norm_files = dir('C:\isbe\dev\background\images\normal512\*.bmp');
bg_mass_idx2 = sort(randsample(1:89, 6));
bg_regions2 = cell(1,6);
for ii = 1:6
    roi = imread(['C:\isbe\dev\background\images\normal1024\', norm_files(bg_mass_idx2(ii)).name]);
    bg_regions2{ii} = roi(1:756, 1:756);
end
save K:\isbe\conferences_and_symposia\iwdm2008\poster\figures\bg_regions2 bg_regions2
%%
%iv) select 12 synthetic masses and shape to be 756x756
syn_mass_idx = randsample(1:101, 12);
syn_files = dir('C:\isbe\dev\mass_model\new_masses\masses080327\*mass*');
%%
for ii = 1:12
    load(['C:\isbe\dev\mass_model\new_masses\masses080327\', syn_files(syn_mass_idx(ii)).name]);
    [r c] = size(mass.subtract_ROI);
    if r > 756
        start_dim = round((r-756) / 2);
        mass.subtract_ROI = mass.subtract_ROI(start_dim:start_dim+755,:);
    end
    if c > 756
        start_dim = round((c-756) / 2);
        mass.subtract_ROI = mass.subtract_ROI(:,start_dim:start_dim+755);
    end
    save(['K:\isbe\conferences_and_symposia\iwdm2008\poster\figures\syn_mass', zerostr(ii,2)], 'mass');
    
end
%%
for ii = 1:6
    imwrite(uint8(bg_regions{ii}), ['K:\isbe\conferences_and_symposia\iwdm2008\poster\figures\bg_mass1_', zerostr(ii,2), '.bmp']);
    %imwrite(bg_regions2{ii}, ['K:\isbe\conferences_and_symposia\iwdm2008\poster\figures\bg_mass2_', zerostr(ii,2), '.bmp']);
end
%%
for ii = 1:6
    write_im_from_colormap(uint8(syn_mass_bg{ii}), ['K:\isbe\conferences_and_symposia\iwdm2008\poster\figures\syn_mass1_', zerostr(ii,2), '.bmp'], gray(256), [10 255]);
    write_im_from_colormap(uint8(syn_norm_bg{ii}), ['K:\isbe\conferences_and_symposia\iwdm2008\poster\figures\syn_mass2_', zerostr(ii,2), '.bmp'], gray(256), [10 255]);
end