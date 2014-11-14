figure; axis equal; hold on;
n = 12;
colors = lines(n);
plot(0,0,'rx');
for ii = 0:n-1
    plot([0 cos(2*ii*pi/n)], [0 sin(2*ii*pi/n)], 'k');
    for jj = 0:n-1
        plot([cos(2*ii*pi/n) cos(2*ii*pi/n)+cos(2*jj*pi/n)], [sin(2*ii*pi/n) sin(2*ii*pi/n)+sin(2*jj*pi/n)], 'color', colors(ii+1,:), 'linestyle', '--');
    end
end
%%
figure; axis equal; hold on;
n = 12;
colors = lines(n);
plot(0,0,'rx');
for ii = 0:n-1
    plot(cos(2*ii*pi/n), sin(2*ii*pi/n), 'ko');
    for jj = 0:n-1
        plot(cos(2*ii*pi/n)+cos(2*jj*pi/n), sin(2*ii*pi/n)+sin(2*jj*pi/n), 'markeredgecolor', colors(ii+1,:), 'marker', 'x');
    end
end
%%
n = 12;
r = (0:10:40)';
c = -(0:10:40)';
rr = [r bsxfun(@plus, r, sin(2*(3:14)*pi/n))];
cc = [c bsxfun(@plus, c, cos(2*(3:14)*pi/n))];

figure; axis equal; hold on; set(gca, 'ydir', 'reverse');
% plot(c, r, 'bo');
% plot(cc, rr, 'rx');
colors = lines(6);
for ii = 1:6
    plot(cc(:,ii+[1 7]), rr(:,ii+[1 7]), 'markeredgecolor', colors(ii,:), 'marker', 'o', 'linestyle', 'none');
end
letters = 'MABCDEFGHIJKLM';
for ii = 1:13
    text(cc(:,ii), rr(:,ii), letters(ii));
end
%%
im = double(imresize(imread([asymmetryroot 'data/contralateral_rfs/misc/o04_010RCC_1024_3797_3365.bmp']), 0.5, 'bilinear'));
n = 12;
num_levels = 4;
% Create DT-CWT of image
dt = dtwavexfm2b(im, num_levels);

r = [256 128]'; %(200:10:240)';
c = [256 300]'; %(200:10:240)';

n_samples = size(r,1);
data = zeros(n_samples, 2*6*num_levels);

rr = [r bsxfun(@plus, r, sin(2*(0.5+(0:11))*pi/n))];
cc = [c bsxfun(@plus, c, cos(2*(0.5+(0:11))*pi/n))];

%Get interpolated dual-tree coefficients
dts = dt_to_pixel_subset(dt, rr, cc);
%
data(:, 1:6*num_levels) = reshape(abs(dt_to_pixel_subset(dt, r, c)), n_samples, []);
%
%for each orientation
for ii = 1:6
    cols = 6*num_levels + 6*(0:num_levels-1) + ii;
    data(:,cols) = angle(squeeze(dts(:,ii+1,1,:) .* conj(dts(:,ii+7,1,:))));

end
%%
num_levels = 4;
for ii = 1:n_samples
    figure;
    for lev = 1:num_levels
        subplot(2,2,lev);
        hold all;
        alim = 0;
        for ori = 1:6
            mag = data(ii,ori+(lev-1)*6);
            phase = data(ii,ori+(num_levels + lev - 1)*6);
            plot([0 mag*cos(phase)], [0 mag*sin(phase)], 'linewidth', 1.5);
            alim = max(alim, mag);
        end
        axis equal;
        %set(gca, 'Xtick', [], 'Ytick', []);
        axis([-alim alim -alim alim]);
        plot([-alim alim], [0 0], 'k--');
        plot([0 0], [-alim alim], 'k--');

    end
end
%%
%Ok dude - let's test this on some synthetic images....
for theta = 15:15:180
    im = create_rect_bar(1, 10, theta, 128, 128, 64, 64);
    
    dt = dtwavexfm2b(im, num_levels);

    r = 64;
    c = 64;

    n_samples = size(r,1);
    data = zeros(n_samples, 2*6*num_levels);

    rr = [r bsxfun(@plus, r, sin(2*(0.5+(0:11))*pi/12))];
    cc = [c bsxfun(@plus, c, cos(2*(0.5+(0:11))*pi/12))];

    %Get interpolated dual-tree coefficients
    dts = dt_to_pixel_subset(dt, rr, cc);
    %
    data(:, 1:6*num_levels) = reshape(abs(dts(:,1,:,:)), n_samples, []);
    %
    %for each orientation
    for ii = 1:6
        cols = 6*num_levels + 6*(0:num_levels-1) + ii;
        data(:,cols) = angle(squeeze(dts(:,ii+1,1,:) .* conj(dts(:,ii+7,1,:))));
    end
    
    ii = 1;
    figure('Name', ['Line angle = ' num2str(theta)]);
    for lev = 1:num_levels
        subplot(2,2,lev);
        hold all;
        alim = 0;
        for ori = 1:6
            mag = data(ii,ori+(lev-1)*6);
            phase = data(ii,ori+(num_levels + lev - 1)*6);
            plot([0 mag*cos(phase)], [0 mag*sin(phase)], 'linewidth', 1.5);
            alim = max(alim, mag);
        end
        legend({'1', '2', '3', '4', '5', '6'}, 'location', 'eastoutside');
        axis equal;
        %set(gca, 'Xtick', [], 'Ytick', []);
        axis([-alim alim -alim alim]);
        plot([-alim alim], [0 0], 'k--');
        plot([0 0], [-alim alim], 'k--');
    end
end
%%
num_levels = 4;
for theta = 61:1:90
    im = create_rect_bar(1, 10, theta, 128, 128, 64, 64);
    
    dt = dtwavexfm2b(im, num_levels+1);

    r = 64;
    c = 64;

    n_samples = size(r,1);

    %Get interpolated dual-tree coefficients
    dts = reshape(dt_to_pixel_subset(dt, r, r), n_samples, 6, num_levels+1);
    mags = reshape(abs(dts(:,:,1:num_levels)), n_samples, []);
    phases = dts(:,:,1:num_levels) .* conj(dts(:,:,2:num_levels+1).^2);
    phases = reshape(atan2(imag(phases), abs(real(phases))), n_samples, []);
    %clear dts;
    
    [dummy max_ori] = max(mags, [], 2);
    max_ori = rem(max_ori-1,6)+1;
    
    %----------------------------------------------------------------------
    for ori = 1:6
        shift_idx = max_ori == ori;
        for lev = 1:num_levels
            cols = 6*(lev-1)+(1:6);
            mags(shift_idx, cols) = circshift(mags(shift_idx, cols), [0 1-ori]);
            phases(shift_idx, cols) = circshift(phases(shift_idx, cols), [0 1-ori]);
        end
    end
    data = [mags phases];
    %clear mags phases;
    %----------------------------------------------------------------------

    ii = 1;
    figure('Name', ['Line angle = ' num2str(theta)]);
    for lev = 1:num_levels
        subplot(2,2,lev);
        hold all;
        alim = 0;
        for ori = 1:6
            mag = data(ii,ori+(lev-1)*6);
            phase = data(ii,ori+(num_levels + lev - 1)*6);
            %plot([0 mag*cos(phase)], [0 mag*sin(phase)], 'linewidth', 1.5);
            plot([ori ori], [0 mag], 'linewidth', 1.5);
            alim = max(alim, mag);
        end
        legend({'1', '2', '3', '4', '5', '6'}, 'location', 'eastoutside');
        axis([0 7 0 alim]);
        %axis equal;
        %axis([-alim alim -alim alim]);
        %plot([-alim alim], [0 0], 'k--');
        %plot([0 0], [-alim alim], 'k--');
    end
end
%%
num_levels = 4;
for theta = 6:17:363
    im = create_gauss_bar(1, 10, theta, 128, 128, 64, 64);
    %im = max(create_rect_bar(1, 10, theta, 128, 128, 64, 64), create_rect_bar(1, 10, theta+90, 128, 128, 64, 64));
    
    dt = dtwavexfm2b(im, num_levels+1);

    r = 64;
    c = 64;

    n_samples = size(r,1);

    %Get interpolated dual-tree coefficients
    dts = reshape(dt_to_pixel_subset(dt, r, r), n_samples, 6, num_levels+1);
    mags = reshape(abs(dts(:,:,1:num_levels)), n_samples, []);
    phases = dts(:,:,1:num_levels) .* conj(dts(:,:,2:num_levels+1).^2);
    phases = reshape(atan2(imag(phases), abs(real(phases))), n_samples, []);
    %clear dts;
    comp = mags.*exp(i*phases);
    
    [max_mags max_idx] = max(mags, [], 2);
    max_ori = rem(max_idx-1,6)+1;
    max_lev = ceil(max_idx/6);
    
    new_mags = zeros(size(mags));
    new_comp = zeros(size(comp));
    %----------------------------------------------------------------------
    for lev = 1:num_levels
        shift_idx = max_lev == lev;
        if any(shift_idx)
            cols = (1:6)+6*(lev-1);
            weights = bsxfun(@rdivide, mags(shift_idx, cols), sum(mags(shift_idx, cols),2));
        
            for lev2 = 1:num_levels
                for ori = 1:6
                    col = 6*(lev2-1)+ori;
                    new_mags(shift_idx, col) = diag(weights * mags(shift_idx, cols).');
                    new_comp(shift_idx, col) = diag(weights * comp(shift_idx, cols).');
                    weights = circshift(weights,[0 1]);
                end
            end
        end
    end

    %data = [new_mags phases];
    data = [new_mags angle(new_comp)];
    %data = [mags phases];
    %clear mags phases;
    %----------------------------------------------------------------------

    ii = 1;
    figure('Name', ['Line angle = ' num2str(theta)]);
    for lev = [2 3]%1:num_levels
        subplot(2,2,lev-1);
        hold all;
        alim = 0;
        for ori = 1:6
            mag = data(ii,ori+(lev-1)*6);
            phase = data(ii,ori+(num_levels + lev - 1)*6);
            plot([0 mag*cos(phase)], [0 mag*sin(phase)], 'linewidth', 1.5);
            alim = max(alim, mag);
        end
        legend({'1', '2', '3', '4', '5', '6'}, 'location', 'eastoutside');
        axis equal;
        axis([-alim alim -alim alim]);
        plot([-alim alim], [0 0], 'k--');
        plot([0 0], [-alim alim], 'k--');
        
        subplot(2,2,lev+1);
        hold all;
        for ori = 1:6
            mag = data(ii,ori+(lev-1)*6);
            phase = data(ii,ori+(num_levels + lev - 1)*6);
            plot([ori ori], [0 mag], 'linewidth', 1.5);
            alim = max(alim, mag);
        end
        legend({'1', '2', '3', '4', '5', '6'}, 'location', 'eastoutside');
        axis([0 7 0 alim]);
    end
end
%%
%im = u_load('C:\isbe\dev\background\images\normal_smooth128\bg001.mat')+create_rect_bar(1, 30, 135, 128, 128, 64, 64);
im = create_gauss_bar(2, 30, 135, 128, 128);
%im = ones(128);
dt = dtwavexfm2(im, 4);
full_dt = dt_to_full_image(dt);
f1 = figure;
f2 = figure;
f3 = figure;
f4 = figure;
for ori = 1:6
    figure(f1);
    subplot(2,3,ori); image(complex2rgb(full_dt(:,:,ori,3))); axis image;
    figure(f2);
    subplot(2,3,ori); imagesc(angle(full_dt(:,:,ori,3))); axis image; colormap(hsv(256));
    figure(f3);
    subplot(2,3,ori); imagesc(angle(full_dt(:,:,ori,3) .* conj(full_dt(:,:,ori,4).^2))); axis image; colormap(hsv(256));
    figure(f4);
    subplot(2,3,ori); imagesc(angle(dt{3}(:,:,ori))); axis image; colormap(hsv(256));
    
end
%%
%im = ones(128);
im = u_load('C:\isbe\dev\background\images\normal_smooth128\bg001.mat')+create_rect_bar(1, 30, 135, 128, 128, 64, 64);
dt = dtwavexfm2(im, 4);
dtb = dtwavexfm2b(im, 4);
full_dt = dt_to_full_image(dt);
full_dtb = dt_to_full_image(dtb);

figure; hold on; 
plot(angle(diag(full_dt(:,:,2,2))), 'r:', 'linewidth', 1.5);
plot(angle(diag(full_dtb(:,:,2,2))), 'g:', 'linewidth', 1.5);

figure; hold on; 
plot(angle(diag(dt{2}(:,:,2))), 'r:', 'linewidth', 1.5);
plot(angle(diag(dtb{2}(:,:,2))), 'g:', 'linewidth', 1.5);
%%
figure;
for ori = 1:6
    subplot(2,3,ori); surf(angle(dt{3}(:,:,ori))); axis image; colormap(hsv(256));    
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%%
centres = (112:143)/2;
for ori = 1:6
    theta = 30*ori - 15;
    
    mag_map = zeros(32);
    phase_map = zeros(32);
    for x = 1:32
        for y = 1:32
            im = create_gauss_blob(2, 1, 128, 128, centres(x), centres(y));
            
            dt = dtwavexfm2b(im, 3);
            mag_map(y, x) = abs(dt{3}(8,8,ori));
            phase_map(y, x) = angle(dt{3}(8,8,ori));
        end
    end
    
    [dummy max_idx] = max(mag_map(:));
    [my mx] = ind2sub([32 32], max_idx);
    
    x_diff = diff(phase_map(1:17,1:17), 1, 2);
    y_diff = diff(phase_map(1:17,1:17), 1, 1);
    
    x_diff(x_diff < -pi) = x_diff(x_diff < -pi) + 2*pi;
    x_diff(x_diff > pi) = x_diff(x_diff > pi) - 2*pi;
    
    y_diff(y_diff < -pi) = y_diff(y_diff < -pi) + 2*pi;
    y_diff(y_diff > pi) = y_diff(y_diff > pi) - 2*pi;
    
    display(['Band ' num2str(ori) ': Maximum located at [' num2str(centres(mx),3) ',' num2str(centres(my),3) '], phase at max location = ', num2str(phase_map(my,mx), 3)]);
    
    figure; imagesc(mag_map); axis image; colormap(jet(256));
    title(['Magnitude map, band ' num2str(ori)]);
    figure; imagesc(phase_map); axis image; colormap(hsv(256)); colorbar;
    title(['Phase map, band ' num2str(ori)]);
    figure; hold all;
    subplot(1,2,1);
    plot(repmat(centres(1:16), 17, 1)', x_diff'); title('X phase gradient');
    subplot(1,2,2);
    plot(repmat(centres(1:16), 17, 1)', y_diff); title('Y phase gradient');
    
end
%%
phases = zeros(128,6);
mags = zeros(128,6);
for x = 1:128
    im = create_gauss_blob(2, 1, 128, 128, 129-x, 60.5);
    dt = dtwavexfm2b(im, 4);
    phases(x,:) = squeeze(angle(dt{3}(8,8,:)))';
    mags(x,:) = squeeze(abs(dt{3}(8,8,:)))';
    clear im dt;
end
im = create_gauss_blob(2, 1, 128, 128, 67.5, 64);
dt = dtwavexfm2b(im, 4);
full_dt = dt_to_full_image(dt, 1:3, 'cubic');
%
f1 = figure;
f2 = figure;
f3 = figure;
f4 = figure;

for ori = 1:6
    figure(f1);
    subplot(2,3,ori); 
    plot(phases(:,ori), 'b'); hold on;
    plot(angle(full_dt(64,:,ori,3)), 'r:'); title(['Band = ' num2str(30*ori - 15) '^{\circ}']);

    figure(f2);
    subplot(2,3,ori); 
    plot(mags(:,ori), 'b'); hold on;
    plot(abs(full_dt(64,:,ori,3)), 'r:'); title(['Band = ' num2str(30*ori - 15) '^{\circ}']);

    figure(f3);
    subplot(2,3,ori); imagesc(angle(full_dt(:,:,ori,3))); axis image; colormap(hsv(256));
    
    figure(f4);
    subplot(2,3,ori); imagesc(abs(full_dt(:,:,ori,3))); axis image; colormap(gray(256));
end
%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
for theta = 7.5:7.5:15
    
    phases = zeros(128,6);
    mags = zeros(128,6);
    for x = 1:128
        im = create_gauss_bar(2, 1, theta, 128, 128, 129-x, 60.5);
        dt = dtwavexfm2b(im, 4);
        phases(x,:) = squeeze(angle(dt{3}(8,8,:)))';
        mags(x,:) = squeeze(abs(dt{3}(8,8,:)))';
        clear im dt;
    end
    im = create_gauss_bar(2, 1, theta, 128, 128, 67.5, 64);
    dt = dtwavexfm2b(im, 4);
    full_dt = dt_to_full_image(dt, 1:3, 'cubic');
    %
    f1 = figure('Name', ['Theta = ', num2str(theta)]);
    f2 = figure('Name', ['Theta = ', num2str(theta)]);
    %f3 = figure('Name', ['Theta = ', num2str(theta)]);
 
    for ori = 1:6
        figure(f1);
        subplot(2,3,ori); 
        plot(phases(:,ori), 'b'); hold on;
        plot(angle(full_dt(64,:,ori,3)), 'r:'); title(['Theta = ' num2str(theta) '^{\circ} Band = ' num2str(30*ori - 15) '^{\circ}']);
        
        figure(f2);
        subplot(2,3,ori); 
        plot(mags(:,ori), 'b'); hold on;
        plot(abs(full_dt(64,:,ori,3)), 'r:'); title(['Theta = ' num2str(theta) '^{\circ} Band = ' num2str(30*ori - 15) '^{\circ}']);
        
        %figure(f3);
        %subplot(2,3,ori); imagesc(angle(full_dt(:,:,ori,3))); axis image; colormap(hsv(256));
    end
    clear im dt full_dt;
end
