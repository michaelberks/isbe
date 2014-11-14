normal_dt = u_load('C:\isbe\dev\background\dual_tree\normal512\o04_008LCC_dual_tree.mat');
mass_dt = u_load('C:\isbe\dev\background\dual_tree\mass_2\mass015_dual_tree.mat');

normal_im = dtwaveifm2(normal_dt, 'near_sym_b','qshift_b');
mass_im = dtwaveifm2(mass_dt, 'near_sym_b','qshift_b');

[rows cols] = size(mass_im);
normal_im_cut = normal_im(1:size(mass_im,1), 1:size(mass_im,2));
normal_dt_cut = dtwavexfm2(normal_im_cut, 5, 'near_sym_b','qshift_b');

figure; imagesc(mass_im); axis image; colormap(gray(256));
figure; imagesc(normal_im_cut); axis image; colormap(gray(256));
%%
normal_dt_cut2 = normal_dt_cut;
mass_dt2 = mass_dt;
for lev = 1:5
    
    normal_dt_cut2{lev} = abs(normal_dt_cut{lev}) .* exp(i*angle(mass_dt{lev}));
    mass_dt2{lev} = abs(mass_dt{lev}) .* exp(i*angle(normal_dt_cut{lev}));
    
end
%
normal_im_cut2 = dtwaveifm2(normal_dt_cut2, 'near_sym_b','qshift_b');
mass_im2 = dtwaveifm2(mass_dt2, 'near_sym_b','qshift_b');

figure; imagesc(mass_im2); axis image; colormap(gray(256));
figure; imagesc(normal_im_cut2); axis image; colormap(gray(256));

%%
lims_1 = [45, 45+rows-1, 57, 57+cols-1];
for lev = 1:4
    lims(lev,:) = ceil(lims_1 / 2^lev); %#ok
end

normal_big_dt = normal_dt;
for lev = 1:4
    normal_big_dt{lev}(lims(lev,1):lims(lev,2), lims(lev,3):lims(lev,4),:) = mass_dt{lev};
end

normal_big_im = dtwaveifm2(normal_big_dt, 'near_sym_b','qshift_b');
figure; imagesc(normal_big_im); axis image; colormap(gray(256));

%%
mammo = u_load('C:\isbe\mammograms\new_CAD\BMP_2004_half\o04_008LCC.mat');
mammo_dt = dtwavexfm2(mammo, 5, 'near_sym_b','qshift_b');
mass_dt = u_load('C:\isbe\dev\background\dual_tree\mass_2\mass015_dual_tree.mat');

lims_1 = [1501, 1501+420-1, 1001, 1001+406-1];
for lev = 1:4
    lims(lev,:) = ceil(lims_1 / 2^lev); %#ok
end

mammo_mass_dt = mammo_dt;
for lev = 1:4
    mammo_mass_dt{lev}(lims(lev,1):lims(lev,2), lims(lev,3):lims(lev,4),:) = mass_dt{lev};
end

mammo_mass_im = dtwaveifm2(mammo_mass_dt, 'near_sym_b','qshift_b');
figure; imagesc(mammo_mass_im); axis image; colormap(gray(256));
