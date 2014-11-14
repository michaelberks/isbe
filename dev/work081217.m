% Take principle components of dual-tree transform co-efficients to build
% models of their distribution

P_dt = cell(4,6);
L_dt = cell(4,6);
B_dt = cell(4,6);
M_dt = cell(4,6);
C_dt = cell(4,6);
mu_dt = cell(4,6);

angle_wrap = [...
    -pi -pi/4 -pi -pi -3*pi/4 -pi/2 ;...
    -pi -pi/4 -pi/4 -pi -3*pi/4 -3*pi/4 ; ...
    -pi -pi/2 -pi/4 -pi -pi -3*pi/4 ;
    -pi -pi/4 -pi -pi -3*pi/4 -pi/2 ;...
    -pi -pi/4 -pi/4 -pi -3*pi/4 -3*pi/4 ; ...
    -pi -pi/2 -pi/4 -pi -pi -3*pi/4 ;...
    ];

for lev = 2:2
    
    %Get a load of data from our texture regions
    [dt_data ilp_data icp_data] = ...
        mb_get_dual_tree_data('C:\isbe\dev\background\dual_tree\normal512', lev, 32, 1, 1);
    
    %Find out which subband is maximal
    [max_dt dt_idx] = max(dt_data, [], 2); clear max_dt
    
    %Circular shift each row of data so maximal sub-band is in column 1
    for sb = 1:1
        band_dt_data = circshift(dt_data(dt_idx == sb, :), [0 1 - sb]);
        band_ilp_data = circshift(ilp_data(dt_idx == sb, :), [0 1 - sb]);
        band_icp_data = circshift(icp_data(dt_idx == sb, :), [0 1 - sb]);
        
        
        %Compute magnitude percentages and angle differences
        dt_lev_mag_percent = abs(band_dt_data(:,2:6)) ./ repmat(abs(band_dt_data(:,1)),1,5);
        ilp_lev_angle_diff = angle(band_ilp_data(:,2:6) .* repmat(conj(band_ilp_data(:,1)),1,5));
        X = zeros(size(band_dt_data,1), 10);
        X(:,1:5) = dt_lev_mag_percent; clear dt_lev_mag_percent;
        X(:,6:10) = ilp_lev_angle_diff; clear ilp_lev_angle_diff;
        
        X_icp = zeros(size(band_dt_data,1), 5);
        for band = 2:6
            
            icp_lev_angle_diff = angle(band_icp_data(:,band) .* conj(band_icp_data(:,1)));
            swap_idx = icp_lev_angle_diff < angle_wrap(sb, band);
            icp_lev_angle_diff(swap_idx) = icp_lev_angle_diff(swap_idx) + 2*pi;
            swap_idx = icp_lev_angle_diff > (angle_wrap(sb, band) + pi);
            icp_lev_angle_diff(swap_idx) = icp_lev_angle_diff(swap_idx) - pi;
            
            X_icp(:, band - 1) = icp_lev_angle_diff; clear icp_lev_angle_diff;
        end
        
        %Compute principal componenets of data
        [P_dt{lev, sb}, L_dt{lev, sb}, B_dt{lev, sb}, M_dt{lev, sb}] = st_pca(X, 0.95);
        [P_icp, L_icp, B_icp, M_icp] = st_pca(X_icp, 1.0);
        %clear X
        
        %So we've got PCA distributions for the ILP angle differences and DT
        %magnitudes - now combine the parameters from these with the maximal ILP
        %angles and DT magnitudes and compute the covariance
        
        X_dt = [abs(band_dt_data(:,1)) angle(band_ilp_data(:,1)) B_icp(:,1) B_dt{lev, sb}];
        C_dt{lev, sb} = cov(X_dt);
        mu_dt{lev, sb} = mean(X_dt);
        %clear X_dt
        
    end
    clear dt_idx
end
%clear band sb lev max_dt

%%        
% Try and rebuild some mammo patches
mammo_list = dir('C:\isbe\dev\background\images\normal512\*.bmp');

rn = randsample(1:length(mammo_list),1);
mammo = double(imread(['C:\isbe\dev\background\images\normal512\', mammo_list(rn).name]));

dt_mammo = dtwavexfm2(mammo, 5, 'near_sym_b','qshift_b');
[ilp_mammo icp_mammo] = mb_dual_tree_transform(dt_mammo);
%
ilp_mammo2 = ilp_mammo;
icp_mammo2 = icp_mammo;
%
for lev = 1:4
    lev_dims = size(dt_mammo{lev}(:,:,1));
    
    %Reshape sub-band
    ilp_data = reshape(ilp_mammo{lev}, [], 6);
    icp_data = reshape(icp_mammo{lev}, [], 6);
    
    [max_mag mag_idx] = max(icp_data, [], 2);
    
    for p = 1:prod(lev_dims);
        sb = mag_idx(p);
        bands = [1:sb-1 sb+1:6];
        
        conditions = ...
            [abs(icp_data(p,sb)) angle(ilp_data(p,sb)) repmat(nan,1,size(P_dt{lev,sb},2))];

        [c_mean c_covar] = ...
            condition_gaussian(mu_dt{lev,sb}, C_dt{lev,sb}, conditions);
        
        new_b = sample_from_normal(c_mean, c_covar, 1);
        
        new_x = M_dt{lev,sb} + new_b * P_dt{lev,sb}';
        
%         new_x = [abs(icp_data(p,bands))./abs(icp_data(p,sb))...
%             angle(ilp_data(p,bands) .*  conj(ilp_data(p,sb)))];
        
        new_mag = abs(icp_data(p,sb)) * new_x(1:5);
        new_phase = angle(ilp_data(p,sb)) + new_x(6:10);

        icp_data(p,bands) = new_mag .* exp(i*angle(icp_data(p,bands)));
        ilp_data(p,bands) = abs(ilp_data(p,bands)) .* exp(i*new_phase);
    end
    
    %Reshape back to original form
    ilp_mammo2{lev} = reshape(ilp_data, lev_dims(1), lev_dims(2), 6);
    icp_mammo2{lev} = reshape(icp_data, lev_dims(1), lev_dims(2), 6);
end

dt_mammo2 = mb_dual_tree_transform_i(ilp_mammo2, icp_mammo2);
dt_mammo2{6} = dt_mammo{6};
mammo2 = dtwaveifm2(dt_mammo2, 'near_sym_b','qshift_b');

display(max(abs(mammo(:) - mammo2(:))));
%figure; image(mammo); axis image; colormap(gray(256));
figure; image(mammo2); axis image; colormap(gray(256));        
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Try and rebuild Lenna

lenna = u_load('C:\isbe\matlab_code\trunk\wavelet_dtcwt_toolbox4_1\lenna.mat');
dt_lenna = dtwavexfm2(lenna, 5, 'near_sym_b','qshift_b');
[ilp_lenna icp_lenna] = mb_dual_tree_transform(dt_lenna);

ilp_lenna2 = ilp_lenna;
%
for lev = 1:4
    lev_dims = size(dt_lenna{lev}(:,:,1));
    
    %Reshape sub-band
    ilp_data = reshape(ilp_lenna{lev}, [], 6);
    icp_data = reshape(icp_lenna{lev}, [], 6);
    
    [max_mag mag_idx] = max(icp_data, [], 2);
    
    for p = 1:prod(lev_dims);
        sb = mag_idx(p);
        bands = [1:sb-1 sb+1:6];
        
        conditions = ...
            [abs(icp_data(p,sb)) angle(ilp_data(p,sb)) repmat(nan,1,size(P_dt{lev,sb},2))];

        [c_mean c_covar] = ...
            condition_gaussian(mu_dt{lev,sb}, C_dt{lev,sb}, conditions);

        new_x = M_dt{lev,sb} + c_mean * P_dt{lev,sb}';
        
%         new_x = [abs(icp_data(p,bands))./abs(icp_data(p,sb))...
%             angle(ilp_data(p,bands) .*  conj(ilp_data(p,sb)))];
        
        new_mag = abs(icp_data(p,sb)) * new_x(1:5);
        new_phase = angle(ilp_data(p,sb)) + new_x(6:10);

        icp_data(p,bands) = new_mag .* exp(i*angle(icp_data(p,bands)));
        ilp_data(p,bands) = abs(ilp_data(p,bands)) .* exp(i*new_phase);
    end
    
    %Reshape back to original form
    ilp_lenna2{lev} = reshape(ilp_data, lev_dims(1), lev_dims(2), 6);
    icp_lenna2{lev} = reshape(icp_data, lev_dims(1), lev_dims(2), 6);
end

dt_lenna2 = mb_dual_tree_transform_i(ilp_lenna2, icp_lenna2);
dt_lenna2{6} = dt_lenna{6};
lenna2 = dtwaveifm2(dt_lenna2, 'near_sym_b','qshift_b');

display(max(abs(lenna(:) - lenna2(:))));
figure; imagesc(lenna); axis image; colormap(gray(256));
figure; imagesc(lenna2); axis image; colormap(gray(256));
 