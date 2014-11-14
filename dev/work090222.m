% Jobs to do!!!!!!!!!!!!!!!!!

% Write a transmit dual_tree from region function

% Look again at the dependence of magnitude on phase, and see if given the
% 6 sub-band phase we can sample matching magnitudes to form good texture

for lev = 2:2
    
    %Get a load of data from our texture regions
    [dt_data ilp_data icp_data] = ...
        mb_get_dual_tree_data('C:\isbe\dev\background\dual_tree\normal512', lev, 32, 1, 1);
    
    %Find out which subband is maximal
    [max_dt dt_idx] = max(dt_data, [], 2); clear max_dt
    
    %Circular shift each row of data so maximal sub-band is in column 1
    for sb = 1:1
        band_dt_data = circshift(dt_data(dt_idx == sb, :), [0 1 - sb]);
    end
    
    %take the principal componenets of the phase
    [P_phase, L_phase, B_phase, M_phase] = st_pca(angle(band_dt_data) / (2*pi), 0.95);
    
     %take the principal componenets of the magnitudes
    [P_mag, L_mag, B_mag, M_mag] = st_pca(abs(band_dt_data) / max(abs(band_dt_data(:))), 0.95);
    
    X_dt = [B_phase B_mag];
    C_dt = cov(X_dt);
    mu_dt = mean(X_dt);
    [RHO,PVAL] = corr(X_dt) %#ok
end