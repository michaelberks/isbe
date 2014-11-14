% % pyr_list = dir([mberksroot, 'background/pyramid/normal512/*pyramid*']);
% % 
% % p_sum_local_pyr = zeros(20, 121);
% % map_sum_local_pyr = zeros(20, 121);
% % 
% % for jj = 1:length(pyr_list)
% %     
% %     display(['processing region ', num2str(jj)]);
% %     
% %     pyr = u_load([mberksroot, 'background/pyramid/normal512/', pyr_list(jj).name]);
% %     
% %     [decomp_corr] = mb_compute_decomp_corr(...
% %         pyr, 'pyr', 'ComputeDecompSubbands', 0, 'ComputeDecompLocal', 11);
% % 
% %     p_sum_local_pyr = p_sum_local_pyr + decomp_corr.local.p_val;
% %     map_sum_local_pyr = map_sum_local_pyr + decomp_corr.local.map;
% % end
% % p_sum_local_pyr = p_sum_local_pyr / jj;
% % 
% % save([mberksroot, 'background/syn/dual_tree/corr_pyr_normal512_local11'], 'map_sum_local_pyr', 'p_sum_local_pyr');

mb_build_dual_tree(...
    'ImageDir', [mberksroot 'new_CAD/BMP_2004_half/'],...
    'OutputDir', [mberksroot 'new_CAD/BMP_2004_dt/'],...
    'ImageFormat', 'mat', 'NumLevels', 6);