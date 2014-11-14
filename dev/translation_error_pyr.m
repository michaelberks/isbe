function translation_errors = translation_error_pyr(region_main, no_oris)  

if ~isa(region_main, 'double')
    region_main = double(region_main);
end

%region_a = region_main(1:512, 1:512);
region_a = region_main;

[pyr_a p_sizes_a] = mb_buildSFpyr(region_a, 5, no_oris-1);
pyramid_a = mb_change_pyramid_form(pyr_a, p_sizes_a);

test_region_a = cell(6, no_oris);
translation_errors.subbands = cell(6, no_oris);
translation_errors.recon = zeros(16);

lims = zeros(6, 2);

for level = 2:6
    lims(level, 1) = (20 * 2^(6-level)) + 1;
    lims(level, 2) = 28 * 2^(6-level);
    
    for ori = 1:no_oris
        test_region_a{level, ori} = pyramid_a{level,ori}...
            (lims(level, 1):lims(level, 2),...
            lims(level, 1):lims(level, 2));
        
        translation_errors.subbands{level, ori} = zeros(16);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute translation errors in each sub-band of the decomposition
% Resolution for translation is at scale of coarsest level
for i = 1:16
    for j = 1:16
        
        region_b = region_main(i*16+(1:512), j*16+(1:512));
        [pyr_b p_sizes_b] = mb_buildSFpyr(region_b, 5, no_oris-1);
        pyramid_b = mb_change_pyramid_form(pyr_b, p_sizes_b);
        
        temp_pyr = pyramid_a;
        
        for level = 2:6
            
            row_offset = i * 2^(6-level);
            col_offset = j * 2^(6-level);
            
            for ori = 1:no_oris
                test_region_b = pyramid_b{level, ori}...
                    ((lims(level,1):lims(level,2)) - row_offset,...
                     (lims(level,1):lims(level,2)) - col_offset);
                 
                temp_pyr{level,ori}...
                    (lims(level,1):lims(level,2),...
                     lims(level,1):lims(level,2)) = test_region_b; 
                
                %translation_errors.subbands{level, ori}(i, j) = ...
                %    sqrt(mean((test_region_a{level,ori}(:)-test_region_b(:)).^2));
                translation_errors.subbands{level, ori}(i, j) = ...
                    mean(abs(test_region_a{level,ori}(:)-test_region_b(:))) /...
                    mean(abs(test_region_a{level,ori}(:)));
                
%                 figure;
%                 subplot(2,2,1); imagesc(pyramid_a{level, ori}); axis image; caxis([-6 6]);
%                 hold on;
%                 plot([lims(level,1) lims(level,2) lims(level,2) lims(level,1) lims(level,1)],...
%                     [lims(level,1) lims(level,1) lims(level,2) lims(level,2) lims(level,1)], 'k-');
%                 subplot(2,2,2); imagesc(pyramid_b{level, ori}); axis image; caxis([-6 6]);
%                 hold on;
%                 plot([lims(level,1) lims(level,2) lims(level,2) lims(level,1) lims(level,1)]-col_offset,...
%                     [lims(level,1) lims(level,1) lims(level,2) lims(level,2) lims(level,1)]-row_offset, 'k-');
%                 subplot(2,2,3); imagesc(test_region_a{level,ori}); axis image; caxis([-5 5]);
%                 subplot(2,2,4); imagesc((test_region_b-test_region_a{level,ori})); axis image; colorbar;
            end
        end
        region_recon = mb_reconSFpyr(mb_change_pyramid_form(temp_pyr), p_sizes_a);
        test_recon_a = region_a(lims(2,1):lims(2,2),lims(2,1):lims(2,2));
        test_recon_b = region_recon(lims(2,1):lims(2,2),lims(2,1):lims(2,2));
        translation_errors.recon(i,j) = ...
            mean(abs(test_recon_a(:)-test_recon_b(:)));
        
%         figure;
%         subplot(1,2,1); imagesc(test_recon_a); axis image; colormap(gray(256));
%         hold on;
%         plot([lims(2,1) lims(2,2) lims(2,2) lims(2,1) lims(2,1)],...
%             [lims(2,1) lims(2,1) lims(2,2) lims(2,2) lims(2,1)], 'k-');
%         subplot(1,2,2); imagesc(test_recon_b); axis image; colormap(gray(256));
%         hold on;
%         plot([lims(2,1) lims(2,2) lims(2,2) lims(2,1) lims(2,1)],...
%             [lims(2,1) lims(2,1) lims(2,2) lims(2,2) lims(2,1)], 'k-');
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute translation errors at level of finest scale
for row_offset = 1:16
    for col_offset = 1:16
        
        region_b = region_main(row_offset+(1:512), col_offset+(1:512));
        [pyr_b p_sizes_b] = mb_buildSFpyr(region_b, 5, no_oris-1);
        pyramid_b = mb_change_pyramid_form(pyr_b, p_sizes_b);
        
        for ori = 1:no_oris
            test_region_b = pyramid_b{2, ori}...
                ((lims(2,1):lims(2,2)) - row_offset,...
                 (lims(2,1):lims(2,2)) - col_offset);

            translation_errors.finest_level{1, ori}(row_offset, col_offset) = ...
                mean(abs(test_region_a{2,ori}(:)-test_region_b(:))) /...
                mean(abs(test_region_a{2,ori}(:)));


        end
        
    end
end

                
        