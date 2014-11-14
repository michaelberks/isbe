cols = repmat(1:256, 256, 1);
rows = cols';
for th = 0:15:165
    cx = 130.5;
    cy = 130.5;
    edge = 1 - create_rect_step(0, 1, th, 256, 256, cx, cy);
    [~, band] = min(abs(th - (15:30:165)));  
    dt = dtwavexfm2b(edge, 3);
    for lev = 2:3
%         [full_tree] = dt_to_pixel_subset(dt, rows, cols, lev);
% 
%         figure;
%         for bb = 1:6
%             subplot(2,3,bb); imgray(complex2rgb(full_tree(65:192,65:192,bb)));
%             plot(cx - 65, cy - 65, 'r+', 'markersize', 10);
%         end

    
        dt_coeffs = squeeze(dt_to_pixel_subset(dt, cy, cx, lev));
        fold_idx = imag(dt_coeffs) < 0;
        dt_coeffs(fold_idx) = conj(dt_coeffs(fold_idx));

        mag = max(abs(dt_coeffs));
        phase = angle(dt_coeffs(band));
        %mag = abs(dt{lev}(33,33,band));
        %phase = abs(dt{lev}(33,33,band));
        display(['Edge feature: Theta = ' num2str(th) ' degs, band = ', num2str(band) ', mag = ' num2str(mag) ', phase = ' num2str(phase) ', level = ' num2str(lev)]);
    end
end
display('');
%%
for th = 0:15:165
    
    cx = 130.5;
    cy = 130.5;
    line = 1 - create_gauss_bar(2, 1, th, 256, 256, cx, cy);
    [~, band] = min(abs(th - (15:30:165)));  
    dt = dtwavexfm2b(line, 3);
    for lev = 2:3
%         [full_tree] = dt_to_pixel_subset(dt, rows, cols, lev);
% 
%         figure;
%         for bb = 1:6
%             subplot(2,3,bb); imgray(complex2rgb(full_tree(65:192,65:192,bb)));
%             plot(cx - 65, cy - 65, 'r+', 'markersize', 10);
%         end

    
        dt_coeffs = squeeze(dt_to_pixel_subset(dt, cy, cx, lev));
        fold_idx = imag(dt_coeffs) < 0;
        dt_coeffs(fold_idx) = conj(dt_coeffs(fold_idx));

        mag = max(abs(dt_coeffs));
        phase = angle(dt_coeffs(band));
        %mag = abs(dt{lev}(33,33,band));
        %phase = abs(dt{lev}(33,33,band));
        display(['Line feature: Theta = ' num2str(th) ' degs, band = ', num2str(band) ', mag = ' num2str(mag,2) ', phase = ' num2str(phase,2) ', level = ' num2str(lev)]);
    end
end
display('');
%%
for th = 0:15:165
    
    cx = 130.5;
    cy = 130.5;
    line = create_gauss_bar(2, 1, th, 256, 256, cx, cy);
    [~, band] = min(abs(th - (15:30:165)));  
    dt = dtwavexfm2b(line, 3);
    for lev = 2:3
        [full_tree] = dt_to_pixel_subset(dt, rows, cols, lev);

        figure;
        for bb = 1:6
            subplot(2,3,bb); imgray(complex2rgb(full_tree(97:160,97:160,bb)));
            plot(cx - 97, cy - 97, 'r+', 'markersize', 10);
        end
    end
end
display('');
%%
cols = repmat(1:256, 256, 1);
rows = cols';
for th = 0:15:75
    cx = 130.5;
    cy = 130.5;
    edge = 1 - create_rect_step(0, 1, th, 256, 256, cx, cy);
    [g2d_responses] = compute_gaussian_2nd_derivatives(edge, [1 2 4 8]);
    
    figure;
    for lev = 1:2
        
        subplot(2,3,3*(lev-1) + 1); imgray(g2d_responses(:,:,lev,1));
        subplot(2,3,3*(lev-1) + 2); imgray(g2d_responses(:,:,lev,2));
        subplot(2,3,3*(lev-1) + 3); imgray(g2d_responses(:,:,lev,3));
        %display(['Edge feature: Theta = ' num2str(th) ' degs, band = ', num2str(band) ', mag = ' num2str(mag) ', phase = ' num2str(phase) ', level = ' num2str(lev)]);
    end
    figure;
    for lev = 3:4
        
        subplot(2,3,3*(lev-3) + 1); imgray(g2d_responses(:,:,lev,1));
        subplot(2,3,3*(lev-3) + 2); imgray(g2d_responses(:,:,lev,2));
        subplot(2,3,3*(lev-3) + 3); imgray(g2d_responses(:,:,lev,3));
        %display(['Edge feature: Theta = ' num2str(th) ' degs, band = ', num2str(band) ', mag = ' num2str(mag) ', phase = ' num2str(phase) ', level = ' num2str(lev)]);
    end
end
%%
for th = 45
    cx = 128;
    cy = 128;
    line = 1-create_rect_bar(1, 0.5, th, 256, 256, cx, cy);
    [g2d_responses] = compute_gaussian_2nd_derivatives(line, [1 2 4 8]);
    clims = [min(g2d_responses(:)) max(g2d_responses(:))];
    
    figure;
    a1 = subplot(1,3,1); hold all;
    a2 = subplot(1,3,2); hold all;
    a3 = subplot(1,3,3); hold all;
    for lev = 1:4
        plot(a1, g2d_responses(cy,95:192,lev,1));
        plot(a2, g2d_responses(cy,95:192,lev,2));
        plot(a3, g2d_responses(cy,95:192,lev,3));
    end
    figure; imgray(line);
    plot([0 256], [cy cy], 'r');
    plot(cx, cy, 'gx');
    %linkaxes([a1 a2 a3]);
end
%
for th = 45
    cx = 128;
    cy = 128;
    edge = create_rect_step(0, 1, th, 256, 256, cx, cy);
    [g2d_responses] = compute_gaussian_2nd_derivatives(edge, [1 2 4 8]);
    clims = [min(g2d_responses(:)) max(g2d_responses(:))];
    
    figure;
    a1 = subplot(1,3,1); hold all;
    a2 = subplot(1,3,2); hold all;
    a3 = subplot(1,3,3); hold all;
    for lev = 1:4
        plot(a1, g2d_responses(cy,95:192,lev,1));
        plot(a2, g2d_responses(cy,95:192,lev,2));
        plot(a3, g2d_responses(cy,95:192,lev,3));
    end
    figure; imgray(edge);
    plot([0 256], [cy cy], 'r');
    plot(cx, cy, 'gx');
    %linkaxes([a1 a2 a3]);
end
%display('');
%%
edge = create_rect_step(0, 1, 45, 256, 256, 128, 128);
[g2d_responses_edge] = compute_gaussian_2nd_derivatives(edge, [1 2 4 8]);
g2d_coeffs_edge = squeeze(g2d_responses_edge(128,129,:,:));
display(g2d_coeffs_edge);

line = 1-create_rect_bar(1, 0.5, 45, 256, 256, 128, 128);
[g2d_responses_line] = compute_gaussian_2nd_derivatives(line, [1 2 4 8]);
g2d_coeffs_line = squeeze(g2d_responses_line(128,128,:,:));
display(g2d_coeffs_line);
