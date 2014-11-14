%**************************************************************************
%**************************************************************************
% Why won't rotation invariance work??
%**************************************************************************

rfs = {'191616', '191681', '191684'};

for ii = 1:3
    num_nodes = zeros(100,1);
    depths = zeros(100,1);
    for jj = 1:100
        tree = u_load(['Z:\asymmetry_project\data\line_detection_rfs\' rfs{ii} '\rf_tree' zerostr(jj,4) '.mat']);
        num_nodes(jj) = length(tree.node);
        depths(jj) = get_tree_depth(tree);
    end
    
    display(['Mean nodes per tree for forest ' rfs{ii} ' = ' num2str(mean(num_nodes))]);
    display(['Mean depth of tree for forest ' rfs{ii} ' = ' num2str(mean(depths))]);
    display(['Max tree depth for forest ' rfs{ii} ' = ' num2str(max(depths))]);
    display('');
end
%%
%Let's create a synthetic CLS, and have a look at the properties of it's
%DT-CWT coefficients as it is rotated about a circle
load C:\isbe\asymmetry_project\data\misc\A.mat
bg = u_load('C:\isbe\asymmetry_project\data\synthetic_backgrounds\smooth512\train\bg00500.mat');
cls = create_gauss_bar(3, 10, 0, 512, 512, 256.5, 256.5);
cls = cls+bg;
phases = zeros(360, 6);
mags = zeros(360, 6);
comps = zeros(360, 6);
mps = zeros(360,12);

mags_aligned = zeros(360, 6);
phases_aligned = zeros(360, 6);
comps_aligned = zeros(360, 6);
mps_aligned = zeros(360, 12);

mags_interp = zeros(360, 6);
phases_interp = zeros(360, 6);
phases_ilp = zeros(360, 6);

lev = 3;
for ii = 1:360
    
        %angled_cls = imrotate(cls, ii,'crop', 'bicubic');
        %angled_cls = angled_cls(129:384, 129:384);
        angled_cls = create_gauss_bar(2, 10, ii, 256, 256, 128.5, 128.5);
        dt = dtwavexfm2b(angled_cls, 4);
        
        dts = dt_to_pixel_subset(dt, 128.5, 128.5);
        
        ilp = dts(:,:,:,lev).*conj(dts(:,:,:,lev+1).^2);
        fold_idx = imag(ilp) < 0;
        ilp(fold_idx) = conj(ilp(fold_idx));
        phases_ilp(ii,:) = squeeze( angle(ilp) )';
        
        fold_idx = imag(dts) < 0;
        dts(fold_idx) = conj(dts(fold_idx));
        
        phases(ii,:) = squeeze(angle(dts(:,:,:,lev)))';
        mags(ii,:) = squeeze(abs(dts(:,:,:,lev)))';
        comps(ii,:) = squeeze(dts(:,:,:,lev)).';
        mps(ii,:) = [mags(ii,:) phases(ii,:)];        
        
        mags_interp(ii,:) = interp1(0:30:180, [mags(ii,6) mags(ii,:)], mod((30:30:180)+ii, 180), 'linear');
        phases_interp(ii,:) = interp1(0:30:180, [phases(ii,6) phases(ii,:)], mod((30:30:180)+ii, 180), 'linear');
        
        mags_aligned(ii,:) = (A_mags(:,:,ii)*mags(ii,:)')';
        phases_aligned(ii,:) = (A_phases(:,:,ii)*phases(ii,:)')';
        comps_aligned(ii,:) = (A_comps(:,:,ii)*comps(ii,:).').';
        mps_aligned(ii,:) = (A_mps(:,:,ii)*mps(ii,:)')';
end
[max_comps idx] = max(comps, [], 2);
max_phases = angle(max_comps);
max_mags = abs(max_comps);

comps_aligned(imag(comps_aligned) < 0) = conj(comps_aligned(imag(comps_aligned) < 0));
% figure; plot(1:360, max_mags);
% axis([1 360 0 max(max_mags(:))]);
% figure; plot(1:360, max_phases);
% axis([1 360 -pi pi]);
%
%%
figure('name', 'DT magnitudes at centre of rotation');
for ii = 1:6
    subplot(2,3,ii); hold all;
    %plot(1:360, mags(:,ii)); axis([1 360 0 max(mags(:))]);
    plot(1:360, mags_interp(:,ii)); axis([1 360 0 max(mags(:))]);
    plot(1:360, mags_aligned(:,ii)); axis([1 360 0 max(mags(:))]);
    plot(1:360, abs(comps_aligned(:,ii))); axis([1 360 0 max(mags(:))]);
    plot(1:360, mps_aligned(:,ii)); axis([1 360 0 max(mags(:))]);
    
    xlabel('Angle (degrees)');
    ylabel('DT-CWT magnitude');
    title(['Band ' num2str(ii)]);
    
    if ii == 1 || ii == 6
        loc = 'southeast';
    else
        loc = 'northeast';
    end
    legend({...
        ['Interpolated, std = ' num2str(std(mags_interp(:,ii)))],...
        ['Magnitude aligned, std = ' num2str(std(mags_aligned(:,ii)))],...
        ['Complex aligned, std = ' num2str(std(abs(comps_aligned(:,ii))))],...
        ['Phase and mag aligned, std = ', num2str(std(mps_aligned(:,ii)))]},...
        'location' , loc);
end

figure('name', 'DT phases at centre of rotation');
for ii = 1:6
    subplot(2,3,ii); hold all;
    plot(1:360, phases_interp(:,ii)); axis([1 360 -pi pi]);
    plot(1:360, phases_aligned(:,ii)); axis([1 360 -pi pi]);
    plot(1:360, angle(comps_aligned(:,ii))); axis([1 360 -pi pi]);
    plot(1:360, mps_aligned(:,ii+6)); axis([1 360 -pi pi]);
    
    xlabel('Angle (degrees)');
    ylabel('DT-CWT phase');
    title(['Band ' num2str(ii)]);
    
    legend({...
        ['Interpolated, std = ' num2str(std(phases_interp(:,ii)))],...
        ['Phase aligned, std = ' num2str(std(phases_aligned(:,ii)))],...
        ['Complex aligned, std = ' num2str(std(angle(comps_aligned(:,ii))))],...
        ['Phase and mag aligned, std = ', num2str(std(mps_aligned(:,ii+6)))]},...
        'location' , 'southeast');
end
%
sort_comps = comps;
sort_ilp = phases_ilp;
for bb = 1:6
    sort_comps(idx == bb,:) = circshift(sort_comps(idx == bb,:), [0 1-bb]);
    sort_ilp(idx == bb,:) = circshift(sort_ilp(idx == bb,:), [0 1-bb]);
end
%
figure('name', 'DT magnitudes at centre of rotation');
for ii = 1:6
    subplot(2,3,ii);
    plot(1:360, abs(sort_comps(:,ii))); axis([1 360 0 max(mags(:))]);
    xlabel('Angle (degrees)');
    ylabel('DT-CWT magnitude');
    title(['Band ' num2str(ii)]);  
end

figure('name', 'DT phases at centre of rotation');
for ii = 1:6
    subplot(2,3,ii);
    plot(1:360, angle(sort_comps(:,ii))); axis([1 360 -pi pi]);
    xlabel('Angle (degrees)');
    ylabel('DT-CWT phase');
    title(['Band ' num2str(ii)]);
end

figure('name', 'DT phases at centre of rotation');
for ii = 1:6
    subplot(2,3,ii);
    plot(1:360, sort_ilp(:,ii)); axis([1 360 -pi pi]);
    xlabel('Angle (degrees)');
    ylabel('DT-CWT phase');
    title(['Band ' num2str(ii)]);
end
%%

%%
figure('name', 'DT magnitudes at centre of rotation');
for ii = 1:6
    subplot(2,3,ii);
    plot(1:360, mags(:,ii)); axis([1 360 0 max(mags(:))]);
    xlabel('Angle (degrees)');
    ylabel('DT-CWT magnitude');
    title(['Band ' num2str(ii)]);  
end

figure('name', 'DT phases at centre of rotation');
for ii = 1:6
    subplot(2,3,ii);
    plot(1:360, phases(:,ii)); axis([1 360 -pi pi]);
    xlabel('Angle (degrees)');
    ylabel('DT-CWT phase');
    title(['Band ' num2str(ii)]);
end
%%
%%
% figure('name', 'DT magnitudes and phases at centre of rotation');
% for ii = 1:6
%     subplot(2,3,ii);
%     polar(phases(:,ii), mags(:,ii));
%     xlabel('Angle (degrees)');
%     ylabel('DT-CWT magnitude');
%     title(['Band ' num2str(ii)]);
% end
%%
figure('name', 'DT phases at centre of rotation');
plot(1:360, phases, '.');
xlabel('Angle (degrees)');
ylabel('DT-CWT phase');
title(['Band ' num2str(ii)]);
legend;
%%

for ii = 1:6
    figure('name', 'DT phases at centre of rotation');
    plot3(1:360, cos(phases(:,ii)), sin(phases(:,ii)));
    xlabel('Angle (degrees)');
    title(['Band ' num2str(ii)]);
end
%%
for ii = 1:6
    figure('name', 'DT phases at centre of rotation');
    plot3(1:360, mags(:,ii).*cos(phases(:,ii)), mags(:,ii).*sin(phases(:,ii))); hold on;
    plot3([1 360], [0 0], [0 0], 'k');
    xlabel('Angle (degrees)');
    title(['Band ' num2str(ii)]);
end
%%
A_mags = zeros(6,6,360);
A_phases = zeros(6,6,360);
A_comps = zeros(6,6,360);
A_mps = zeros(12,12,360);
for jj = 1:360
    mags_theta = zeros(6,200);
    mags_zero = zeros(6,200);
    
    phases_theta = zeros(6,200);
    phases_zero = zeros(6,200);

    comps_theta = zeros(6,200);
    comps_zero = zeros(6,200);

    mps_theta = zeros(12,200);
    mps_zero = zeros(12,200);

    for ii = 1:200
        bg = u_load(['C:\isbe\asymmetry_project\data\synthetic_backgrounds\smooth512\train\bg' zerostr(ii,5) '.mat']);
        cls = create_gauss_bar(2, 10, 0, 512, 512, 256.5, 256.5);
        cls = cls+bg;

        cls_zero = cls(129:384, 129:384);
        dt_zero = dtwavexfm2b(cls_zero, 4);
        dts_zero = dt_to_pixel_subset(dt_zero, 128.5, 128.5);
        fold_idx = imag(dts_zero) < 0;
        dts_zero(fold_idx) = conj(dts_zero(fold_idx));
        
        mags_zero(:,ii) = squeeze(abs(dts_zero(:,:,:,3)));
        phases_zero(:,ii) = squeeze(angle(dts_zero(:,:,:,3)));
        comps_zero(:,ii) = squeeze(dts_zero(:,:,:,3)).';
        mps_zero(:,ii) = [mags_zero(:,ii); phases_zero(:,ii)];
    
        cls_theta = imrotate(cls, jj, 'crop', 'bicubic');
        cls_theta = cls_theta(129:384, 129:384);
        dt_theta = dtwavexfm2b(cls_theta, 4);
        dts_theta = dt_to_pixel_subset(dt_theta, 128.5, 128.5);
        fold_idx = imag(dts_theta) < 0;
        dts_theta(fold_idx) = conj(dts_theta(fold_idx));
        
        mags_theta(:,ii) = squeeze(abs(dts_theta(:,:,:,3)));
        phases_theta(:,ii) = squeeze(angle(dts_theta(:,:,:,3)));
        comps_theta(:,ii) = squeeze(dts_theta(:,:,:,3)).';
        mps_theta(:,ii) = [mags_theta(:,ii); phases_theta(:,ii)];
    end

    A_mags(:,:,jj) = mags_zero / mags_theta;
    A_phases(:,:,jj) = phases_zero / phases_theta;
    A_comps(:,:,jj) = comps_zero / comps_theta;
    A_mps(:,:,jj) = mps_zero / mps_theta;
end
save C:\isbe\asymmetry_project\data\misc\A.mat A*
%%
for ii =1:6
    figure;
    for jj = 1:6
        subplot(2,3,jj); plot(1:360, squeeze(A_mags(ii,jj,:)));
    end
end
for ii =1:6
    figure;
    for jj = 1:6
        subplot(2,3,jj); plot(1:360, squeeze(A_phases(ii,jj,:)));
    end
end
%%
for ii = 0:12:348
    figure; colormap(jet(256));
    for jj = 1:12
        subplot(3,4,jj);
        imagesc(A(:,:,ii+jj)); axis image; caxis([min(A(:)) max(A(:))]);
        title(['Angle ' num2str(ii+jj)]);
    end
end
%%
mags_theta = zeros(6,200);
mags_zero = zeros(6,200);

phases_theta = zeros(6,200);
phases_zero = zeros(6,200);

comps_theta = zeros(6,200);
comps_zero = zeros(6,200);

mps_theta = zeros(12,200);
mps_zero = zeros(12,200);

for ii = 1:200
    bg = u_load(['C:\isbe\asymmetry_project\data\synthetic_backgrounds\smooth512\train\bg' zerostr(ii,5) '.mat']);
    cls = create_gauss_bar(2, 10, 0, 512, 512, 256.5, 256.5);
    cls = cls+bg;

    cls_zero = cls(129:384, 129:384);
    dt_zero = dtwavexfm2b(cls_zero, 4);
    dts_zero = dt_to_pixel_subset(dt_zero, 128.5, 128.5);

    mags_zero(:,ii) = squeeze(abs(dts_zero(:,:,:,3)));
    phases_zero(:,ii) = squeeze(angle(dts_zero(:,:,:,3)));
    comps_zero(:,ii) = squeeze(dts_zero(:,:,:,3)).';
    mps_zero(:,ii) = [mags_zero(:,ii); phases_zero(:,ii)];
    
    cls_theta = imrotate(cls, 180, 'crop', 'bicubic');
    cls_theta = cls_theta(129:384, 129:384);
    dt_theta = dtwavexfm2b(cls_theta, 4);
    dts_theta = dt_to_pixel_subset(dt_theta, 128.5, 128.5);

    mags_theta(:,ii) = squeeze(abs(dts_theta(:,:,:,3)));
    phases_theta(:,ii) = squeeze(angle(dts_theta(:,:,:,3)));
    comps_theta(:,ii) = squeeze(dts_theta(:,:,:,3)).';
    mps_theta(:,ii) = [mags_theta(:,ii); phases_theta(:,ii)];
end
A_mags = mags_zero / mags_theta;
A_phases = phases_zero / phases_theta;
A_comps = comps_zero / comps_theta;
A_mps = mps_zero / mps_theta;
%%
figure('name', 'DT magnitudes at centre of rotation');
for ii = 1:6
    subplot(2,3,ii);
    plot(1:200, mags_zero(ii,:));
    xlabel('Angle (degrees)');
    ylabel('DT-CWT magnitude');
    title(['Band ' num2str(ii)]);
end
%%
[roc_mag0, auc_mag0] = compute_roc_image_set_lines(...
    'C:\isbe\asymmetry_project\data\synthetic_lines\lines512\',...
    'C:\isbe\asymmetry_project\data\synthetic_lines\lines512\results\191819\', 'centre_line');
[roc_phase0, auc_phase0] = compute_roc_image_set_lines(...
    'C:\isbe\asymmetry_project\data\synthetic_lines\lines512\',...
    'C:\isbe\asymmetry_project\data\synthetic_lines\lines512\results\191820\', 'centre_line');
[roc_phase1, auc_phase1] = compute_roc_image_set_lines(...
    'C:\isbe\asymmetry_project\data\synthetic_lines\lines512\',...
    'C:\isbe\asymmetry_project\data\synthetic_lines\lines512\results\191821\', 'centre_line');
[roc_mag1, auc_mag1] = compute_roc_image_set_lines(...
    'C:\isbe\asymmetry_project\data\synthetic_lines\lines512\',...
    'C:\isbe\asymmetry_project\data\synthetic_lines\lines512\results\191822\', 'centre_line');
%
[roc_all0, auc_all0] = compute_roc_image_set_lines(...
    'C:\isbe\asymmetry_project\data\synthetic_lines\lines512\',...
    'C:\isbe\asymmetry_project\data\synthetic_lines\lines512\results\191838\', 'centre_line');
[roc_all1, auc_all1] = compute_roc_image_set_lines(...
    'C:\isbe\asymmetry_project\data\synthetic_lines\lines512\',...
    'C:\isbe\asymmetry_project\data\synthetic_lines\lines512\results\191839\', 'centre_line');

[roc_clock0, auc_clock0] = compute_roc_image_set_lines(...
    'C:\isbe\asymmetry_project\data\synthetic_lines\lines512\',...
    'C:\isbe\asymmetry_project\data\synthetic_lines\lines512\results\191681\', 'centre_line');
[roc_clock1, auc_clock1] = compute_roc_image_set_lines(...
    'C:\isbe\asymmetry_project\data\synthetic_lines\lines512\',...
    'C:\isbe\asymmetry_project\data\synthetic_lines\lines512\results\191844\', 'centre_line');
%


[roc_all3, auc_all3] = compute_roc_image_set_lines(...
    'C:\isbe\asymmetry_project\data\synthetic_lines\lines512\',...
    'C:\isbe\asymmetry_project\data\synthetic_lines\lines512\results\191616\', 'centre_line');
%%
[roc_conj1, auc_conj1] = compute_roc_image_set_lines(...
    'C:\isbe\asymmetry_project\data\synthetic_lines\lines512\',...
    'C:\isbe\asymmetry_project\data\synthetic_lines\lines512\results\191865\', 'centre_line');

%%
[roc_conj, auc_conj] = compute_roc_image_set_lines(...
    'C:\isbe\asymmetry_project\data\synthetic_lines\old_lines\lines512\',...
    'C:\isbe\asymmetry_project\data\synthetic_lines\old_lines\lines512\results\191856\', 'centre_line');
[roc_conj3, auc_conj3] = compute_roc_image_set_lines(...
    'C:\isbe\asymmetry_project\data\synthetic_lines\old_lines\lines512\',...
    'C:\isbe\asymmetry_project\data\synthetic_lines\old_lines\lines512\results\191866\', 'centre_line');
%%

figure;
hold all; axis equal; axis([0 1 0 1]); 
plot(roc_all0(:,1), roc_all0(:,2), '--', 'linewidth', 2);
plot(roc_all1(:,1), roc_all1(:,2), '--', 'linewidth', 2);
plot(roc_phase0(:,1), roc_phase0(:,2), '--', 'linewidth', 2);
plot(roc_phase1(:,1), roc_phase1(:,2), '--', 'linewidth', 2);
plot(roc_mag0(:,1), roc_mag0(:,2), '--', 'linewidth', 2);
plot(roc_mag1(:,1), roc_mag1(:,2), '--', 'linewidth', 2);
plot(roc_clock0(:,1), roc_clock0(:,2), '--', 'linewidth', 2);
plot(roc_clock1(:,1), roc_clock1(:,2), '--', 'linewidth', 2);
plot(roc_conj(:,1), roc_conj(:,2), '--', 'linewidth', 2);
plot(roc_all3(:,1), roc_all3(:,2), '--', 'linewidth', 2);
legend({...
    ['All AUC = ' num2str(auc_all0)], ['All rotated AUC = ' num2str(auc_all1)],...
    ['Phase AUC = ' num2str(auc_phase0)], ['Phase rotated AUC = ' num2str(auc_phase1)],...
    ['Mag AUC = ' num2str(auc_mag0)], ['Mag rotated AUC = ' num2str(auc_mag1)],...
    ['Clock AUC = ' num2str(auc_clock0)], ['Clock rotated AUC = ' num2str(auc_clock1)],...
    ['Conj AUC = ' num2str(auc_conj)], ['All 3x3 AUC = ' num2str(auc_all3)]});
%%
interp_coeffs = zeros(180, 6);
for ii = 1:180
    interp_coeffs(ii,:) = interp1(0:30:180, [0 1 0 0 0 0 0 ], mod((0:30:150)+ii, 180), 'linear');
end