function translation_errors = tranlsation_error_dual_tree(region_main)  

if ~isa(region_main, 'double')
    region_main = double(region_main);
end

%Take sample region from main image to act as source
region_a = region_main(1:512, 1:512);
%region_a = region_main;

%Compute Dual-Tree decomposition of source region
[dtl_a, dth_a] = dtwavexfm2(region_a,5,'near_sym_b','qshift_b');

%Extrac test regions from the decomposition and save limits
test_region_a = cell(5, 6);
translation_errors.subbands = cell(1, 6);
translation_errors.recon = zeros(8);
lims = zeros(5, 2);

for level = 1:5
    %Limits approproately downscaled to account for sub-sampling
    lims(level, 1) = (20 * 2^(4-level)) + 1;
    lims(level, 2) = 28 * 2^(4-level);
    
    for ori = 1:6
        test_region_a{level, ori} = dth_a{level}...
            (lims(level, 1):lims(level, 2),...
             lims(level, 1):lims(level, 2),ori);
         
        %Pre-allocate errory array
        translation_errors.subbands{level, ori} = zeros(8);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute translation errors in each sub-band of the decomposition
% Resolution for translation is at scale of coarsest level
for i = 1:16
    for j = 1:16
        
        region_b = region_main(i*16+(1:512), j*16+(1:512));
        
        [dtl_b, dth_b] = dtwavexfm2(region_b,5,'near_sym_b','qshift_b');
        
        temp_dth = dth_a;
        
        for level = 1:4
            
            row_offset = i * 2^(4-level);
            col_offset = j * 2^(4-level);
            
%             display(['row offset = ', num2str(row_offset)]);
%             display(['row offset = ', num2str(row_offset)]);
%             display(lims(level,:));
            
            for ori = 1:6
                test_region_b = dth_b{level}...
                    ((lims(level,1):lims(level,2)) - row_offset,...
                     (lims(level,1):lims(level,2)) - col_offset, ori);
                 
                 temp_dth{level}...
                    (lims(level,1):lims(level,2),...
                     lims(level,1):lims(level,2), ori) = test_region_b;
                
                %shift_invariant_errors{level, ori}(i, j) = ...
                %    sqrt(mean((test_region_a{level,ori}(:)-test_region_b(:)).^2));
                translation_errors.subbands{level, ori}(i, j) = ...
                    mean(abs(test_region_a{level,ori}(:)-test_region_b(:))) /...
                    mean(abs(test_region_a{level,ori}(:)));
                
%                 figure;
%                 subplot(2,2,1); imagesc(abs(dth_a{level}(:,:,ori))); axis image; colorbar;
%                 hold on;
%                 plot([lims(level,1) lims(level,2) lims(level,2) lims(level,1) lims(level,1)],...
%                     [lims(level,1) lims(level,1) lims(level,2) lims(level,2) lims(level,1)], 'k-');
%                 subplot(2,2,2); imagesc(abs(dth_b{level}(:,:,ori))); axis image; colorbar;
%                 hold on;
%                 plot([lims(level,1) lims(level,2) lims(level,2) lims(level,1) lims(level,1)]-col_offset,...
%                     [lims(level,1) lims(level,1) lims(level,2) lims(level,2) lims(level,1)]-row_offset, 'k-');
%                 subplot(2,2,3); imagesc(abs(test_region_a{level,ori})); axis image; colorbar
%                 subplot(2,2,4); imagesc(abs(test_region_b-test_region_a{level,ori})); axis image; colorbar;
            end
        end
        region_recon = dtwaveifm2(dtl_a, temp_dth,'near_sym_b','qshift_b');
        test_recon_a = region_a(321:448,321:448);
        test_recon_b = region_recon(321:448,321:448);
        translation_errors.recon(i,j) = ...
            mean(abs(test_recon_a(:)-test_recon_b(:)));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute translation errors at level of finest scale
for row_offset = 1:16
    for col_offset = 1:16
        
        region_b = region_main(row_offset+(1:512), col_offset+(1:512));
        
        [dtl_b, dth_b] = dtwavexfm2(region_b,5,'near_sym_b','qshift_b');
        
        for ori = 1:6
            
            test_region_b = dth_b{1}...
                    ((lims(1,1):lims(1,2)) - row_offset,...
                     (lims(1,1):lims(1,2)) - col_offset, ori);
            

            translation_errors.finest_level{1, ori}(row_offset, col_offset) = ...
                mean(abs(test_region_a{1,ori}(:)-test_region_b(:))) /...
                mean(abs(test_region_a{1,ori}(:)));

        end
        
    end
end

                
        