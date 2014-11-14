function display_full_tree_info(image_in, levels)

if nargin < 2
    levels = 6;
end

dual_tree = dtwavexfm2(image_in, levels);
full_tree = dt_to_full_image(dual_tree);
clear dual_tree;
%
f1 = figure; imagesc(image_in); axis image; colormap(gray(256));

colors = lines(levels);
bands = [1 2 3 6 5 4];
while 1
    
    [xi,yi,P] = impixel; %#ok
    
    if ~isempty(xi)

        coeffs = squeeze(full_tree(yi,xi,:,:)); %dim 1 = bands, dim 2 = levels
        
        figure(...
            'Name', ['Interpolated dual-tree coefficients at (x,y) = (', num2str(xi), ',',num2str(yi), ')']);
        for band = 1:6
            subplot(2,3,bands(band));
            compass(max(coeffs(:)), 'w'); hold on;
            %compass(coeffs(band,:)); 
            for level = 1:levels
                h(level) = plot([0 coeffs(band,level)], 'Color', colors(level,:), 'LineWidth', 2); %#ok
                legend_string(level,:) =  ['Level ', num2str(level)]; %#ok
            end
            if band == 1;
                legend(h, legend_string); clear h legend_string;
            end
            title(['Band ', num2str(band)]);
        end
        
        for level = 1:levels-1
            ilp(:,level) = abs(coeffs(:,level)) .* exp(i*(angle(coeffs(:,level)) - 2*angle(coeffs(:,level+1)))); %#ok
        end
        idx = real(ilp) < 0;
        ilp(idx) = complex(-real(ilp(idx)), imag(ilp(idx))); %#ok
        
        figure(...
            'Name', ['ILP coefficients at (x,y) = (', num2str(xi), ',',num2str(yi), ')']);
        for band = 1:6
            subplot(2,3,bands(band));
            compass(max(coeffs(:)), 'w'); hold on;
            %compass(ilp(band,:));
            for level = 1:levels-1
                h(level) = plot([0 ilp(band,level)], 'Color', colors(level,:), 'LineWidth', 2); %#ok
                legend_string(level,:) =  ['Level ', num2str(level)]; %#ok
            end
            if band == 1
                legend(h, legend_string); clear h legend_string;
            end
            title(['Band ', num2str(band)]);
        end
        set(0,'CurrentFigure', f1);
    else
        break;
    end
end   
    
