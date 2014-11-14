%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to test the transform of a dual-tree into ILP and ICP coefficients
% and then reconstruct the original tree

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Let's start with Lenna
lenna = u_load('C:\isbe\matlab_code\trunk\wavelet_dtcwt_toolbox4_1\lenna.mat');

% compute dual-tree for each image
dt_lenna = dtwavexfm2(lenna, 5, 'near_sym_b','qshift_b');

%Transform to ILP and ICP
[ilp_lenna icp_lenna] = mb_dual_tree_transform(dt_lenna);
%
%Transform back again and compare
[dt_lenna_i] = mb_dual_tree_transform_i(ilp_lenna, icp_lenna);
for lev = 1:5
    for band = 1:6
        display(['Level ', num2str(lev),...
            ' Sum of reconstruction erros = ', ...
            num2str(sum(sum(abs(dt_lenna{lev}(:,:,band) - dt_lenna_i{lev}(:,:,band))))),...
            ' Maximum reconstruction error = ', ...
            num2str(max(max(abs(dt_lenna{lev}(:,:,band) - dt_lenna_i{lev}(:,:,band)))))]);
    end
end
%clear
%Bonza, we have perfect reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Now lets check rotational invariance
lenna = u_load('C:\isbe\matlab_code\trunk\wavelet_dtcwt_toolbox4_1\lenna.mat');
lenna90 = imrotate(lenna, 90);
lenna180 = imrotate(lenna, 180);
lenna270 = imrotate(lenna, 270);

% compute dual-tree for each image
dt_lenna = dtwavexfm2(lenna, 5, 'near_sym_b','qshift_b');
dt_lenna90 = dtwavexfm2(lenna90, 5, 'near_sym_b','qshift_b');
dt_lenna180 = dtwavexfm2(lenna180, 5, 'near_sym_b','qshift_b');
dt_lenna270 = dtwavexfm2(lenna270, 5, 'near_sym_b','qshift_b');

%
% compute ILP and ICP coefficients for each image
[ilp_0 icp_0] = mb_dual_tree_transform(dt_lenna);
[ilp_90 icp_90] = mb_dual_tree_transform(dt_lenna90);
[ilp_180 icp_180] = mb_dual_tree_transform(dt_lenna180);
[ilp_270 icp_270] = mb_dual_tree_transform(dt_lenna270);
%%
for lev = 1:4
    for ori = 1:3
        
        ilp_diff_90{lev}(:,:,ori) = ilp_0{lev}(:,:,ori) - rot90(ilp_90{lev}(:,:,ori+3),-1);
        ilp_diff_180{lev}(:,:,ori) = ilp_0{lev}(:,:,ori) - rot90(ilp_180{lev}(:,:,ori),-2);
        ilp_diff_270{lev}(:,:,ori) = ilp_0{lev}(:,:,ori) - rot90(ilp_270{lev}(:,:,ori+3),-3);
        
        icp_diff_90{lev}(:,:,ori) = icp_0{lev}(:,:,ori) -...
            rot90(complex(imag(icp_90{lev}(:,:,ori+3)), -real(icp_90{lev}(:,:,ori+3))),-1);
        icp_diff_180{lev}(:,:,ori) = icp_0{lev}(:,:,ori) -...
            rot90(-icp_180{lev}(:,:,ori),-2);
        icp_diff_270{lev}(:,:,ori) = icp_0{lev}(:,:,ori) -...
            rot90(complex(-imag(icp_270{lev}(:,:,ori+3)), real(icp_270{lev}(:,:,ori+3))),-3);
        
    end
    for ori = 4:6
        
        ilp_diff_90{lev}(:,:,ori) = ilp_0{lev}(:,:,ori) - rot90(ilp_90{lev}(:,:,ori-3),-1);
        ilp_diff_180{lev}(:,:,ori) = ilp_0{lev}(:,:,ori) - rot90(ilp_180{lev}(:,:,ori),-2);
        ilp_diff_270{lev}(:,:,ori) = ilp_0{lev}(:,:,ori) - rot90(ilp_270{lev}(:,:,ori-3),-3);
        
        icp_diff_90{lev}(:,:,ori) = icp_0{lev}(:,:,ori) -...
            rot90(complex(imag(icp_90{lev}(:,:,ori-3)), -real(icp_90{lev}(:,:,ori-3))),-1);
        icp_diff_180{lev}(:,:,ori) = icp_0{lev}(:,:,ori) -...
            rot90(-icp_180{lev}(:,:,ori),-2);
        icp_diff_270{lev}(:,:,ori) = icp_0{lev}(:,:,ori) -...
            rot90(complex(-imag(icp_270{lev}(:,:,ori-3)), real(icp_270{lev}(:,:,ori-3))),-3);
        
    end
    display(['Level ', num2str(lev),...
        ' Sum of 90 degree ilp erros = ', ...
        num2str(sum(abs(ilp_diff_90{lev}(:)))),...
        ' Sum of 90 degree icp erros = ', ...
        num2str(sum(abs(icp_diff_90{lev}(:))))]);
    display(['Level ', num2str(lev),... 
        ' Sum of 180 degree ilp erros = ', ...
        num2str(sum(abs(ilp_diff_180{lev}(:)))),...
        ' Sum of 180 degree icp erros = ', ...
        num2str(sum(abs(icp_diff_180{lev}(:))))]);
    display(['Level ', num2str(lev),...
        ' Sum of 270 degree ilp erros = ', ...
        num2str(sum(abs(ilp_diff_270{lev}(:)))),...
        ' Sum of 270 degree icp erros = ', ...
        num2str(sum(abs(icp_diff_270{lev}(:)))),...
        ]);
    
    max_icp0 = max(icp_0{lev}, [], 3);
    max_icp90 = max(icp_90{lev}, [], 3);
    max_icp180 = max(icp_180{lev}, [], 3);
    max_icp270 = max(icp_270{lev}, [], 3);
    
    max_icp_diff_90{lev} = max_icp0 - rot90(complex(imag(max_icp90), -real(max_icp90)),-1);
    max_icp_diff_180{lev} = max_icp0 - rot90(-max_icp180,-2);
    max_icp_diff_270{lev} = max_icp0 - rot90(complex(-imag(max_icp270), real(max_icp270)),-3);
    
    display(['Level ', num2str(lev),...
        ' Sum of 90 degree max icp erros = ', ...
        num2str(sum(abs(max_icp_diff_90{lev}(:))))]);
    display(['Level ', num2str(lev),... 
        ' Sum of 180 degree max icp erros = ', ...
        num2str(sum(abs(max_icp_diff_180{lev}(:))))]);
    display(['Level ', num2str(lev),...
        ' Sum of 270 degree max icp erros = ', ...
        num2str(sum(abs(max_icp_diff_270{lev}(:)))),...
        ]);

end
%%
% A change=======
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The problem occurs when the maximal band differs for the dual-tree, ilp
% and icp coefficients - lets see how many location this is true for

% Make a structure consisting of the index to the maximal band in each of
% the dt, ilp and icp coeffs (only need do this 0 degree image)
max_idx = cell(4,1);
diff_idx12 = cell(4,1);
diff_idx13 = cell(4,1);
diff_idx23 = cell(4,1);

for lev = 1:4
    [r c] = size(dt_lenna{lev}(:,:,1));
    max_idx{lev} = zeros(r, c, 3);
    [m_dt max_idx{lev}(:,:,1)] = max(dt_lenna{lev}, [], 3);
    [m_ilp max_idx{lev}(:,:,2)] = max(ilp_0{lev}, [], 3);
    [m_icp max_idx{lev}(:,:,3)] = max(icp_0{lev}, [], 3);
    clear dummy r c;
    
    diff_idx12{lev} = abs(max_idx{lev}(:,:,1) - max_idx{lev}(:,:,2));
    diff_idx13{lev} = abs(max_idx{lev}(:,:,1) - max_idx{lev}(:,:,3));
    diff_idx23{lev} = abs(max_idx{lev}(:,:,2) - max_idx{lev}(:,:,3));
    
    diff_idx12{lev}(diff_idx12{lev} == 4) = 2;
    diff_idx12{lev}(diff_idx12{lev} == 5) = 1;
    diff_idx13{lev}(diff_idx13{lev} == 4) = 2;
    diff_idx13{lev}(diff_idx13{lev} == 5) = 1;
    diff_idx23{lev}(diff_idx23{lev} == 4) = 2;
    diff_idx23{lev}(diff_idx23{lev} == 5) = 1;
    
    % Look at the counts for each type of difference
%     [sum(diff_idx12{lev}(:)>0) sum(~diff_idx12{lev}(:));...
%      sum(diff_idx13{lev}(:)>0) sum(~diff_idx13{lev}(:));...
%      sum(diff_idx23{lev}(:)>0) sum(~diff_idx23{lev}(:))]
 
    %Do difference occur in small or large magnitude coeffs
    display([mean(abs(m_dt(diff_idx12{lev}(:)>0))), mean(abs(m_dt(~diff_idx12{lev}(:))))]);
    display([mean(abs(m_dt(diff_idx13{lev}(:)>0))), mean(abs(m_dt(~diff_idx13{lev}(:))))]);
    display([mean(abs(m_dt(diff_idx23{lev}(:)>0))), mean(abs(m_dt(~diff_idx23{lev}(:))))]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Who or what to trust????
%Look at the ICP angles across sub-bands
%%
for lev = 1:4
    figure; hold on; axis equal;
    colors = 'rgbmck';
    for band = 1:6
        quiver(real(icp_0{lev}(:,:,band)), imag(icp_0{lev}(:,:,band)), colors(band));
    end    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Try calculating the angle difference between ICP coefficient in each band
% and the angle of the band - we'll then look at the angle differences for
% each of maximal bands for dt/ilp/icp magnitudes

band_angles = [pi/10 pi/4 2*pi/5 3*pi/5 3*pi/4 9*pi/10];
band_angle_lims = [0, 7*pi/40 13*pi/40 20*pi/40 27*pi/40 33*pi/40 pi]; 
angle_diff = cell(4,1);
angle_map = cell(4,1);
angle_diff_dt = cell(4,1);
angle_diff_ilp = cell(4,1);
angle_diff_icp = cell(4,1);
angle_map_dt = cell(4,1);
angle_map_ilp = cell(4,1);
angle_map_icp = cell(4,1);

for lev = 1:4
    angle_diff{lev} = zeros(size(dt_lenna{lev}));
    
    ad_dt = zeros(size(dt_lenna{lev}(:,:,1)));
    ad_ilp = zeros(size(dt_lenna{lev}(:,:,1)));
    ad_icp = zeros(size(dt_lenna{lev}(:,:,1)));
    am_dt = zeros(size(dt_lenna{lev}(:,:,1)));
    am_ilp = zeros(size(dt_lenna{lev}(:,:,1)));
    am_icp = zeros(size(dt_lenna{lev}(:,:,1)));
    
    for band = 1:6
        %calculate the angle difference
        ad = angle(icp_0{lev}(:,:,band)) - band_angles(band);
        
        %take the difference modulus pi, so angles in range (-pi/2, pi/2)
        ad = mb_mod(ad, pi);
        angle_diff{lev}(:,:,band) = ad;
        
        ad_dt(max_idx{lev}(:,:,1)==band) = ad(max_idx{lev}(:,:,1)==band);
        ad_ilp(max_idx{lev}(:,:,2)==band) = ad(max_idx{lev}(:,:,2)==band);
        ad_icp(max_idx{lev}(:,:,3)==band) = ad(max_idx{lev}(:,:,3)==band);
        
        %get logical map for angles outside range of band
        a_map = band_angle_lims(band) <= angle(icp_0{lev}(:,:,band)) & ...
            angle(icp_0{lev}(:,:,band)) < band_angle_lims(band+1);
        angle_map{lev}(:,:,band) = a_map;
        
        am_dt(max_idx{lev}(:,:,1)==band) = a_map(max_idx{lev}(:,:,1)==band);
        am_ilp(max_idx{lev}(:,:,2)==band) = a_map(max_idx{lev}(:,:,2)==band);
        am_icp(max_idx{lev}(:,:,3)==band) = a_map(max_idx{lev}(:,:,3)==band);
        
    end
    
    angle_diff_dt{lev} = ad_dt;
    angle_diff_ilp{lev} = ad_ilp;
    angle_diff_icp{lev} = ad_icp;
    angle_map_dt{lev} = am_dt;
    angle_map_ilp{lev} = am_ilp;
    angle_map_icp{lev} = am_icp;
    
end
%%
for lev = 1:4
    figure('name', 'Dual-tree'); imagesc(angle_diff_dt{lev}); axis image; colormap(hsv); caxis([-pi/2, pi/2]);
    figure('name', 'ILP'); imagesc(angle_diff_ilp{lev}); axis image; colormap(hsv); caxis([-pi/2, pi/2]);
    figure('name', 'ICP'); imagesc(angle_diff_icp{lev}); axis image; colormap(hsv); caxis([-pi/2, pi/2]);
end
%%
for lev = 1:4
    figure('name', 'Dual-tree'); imagesc(angle_map_dt{lev}); axis image; 
    figure('name', 'ILP'); imagesc(angle_map_ilp{lev}); axis image; 
    figure('name', 'ICP'); imagesc(angle_map_icp{lev}); axis image; 
end
%%
for lev = 1:4
    display(['Number of mis-aligned points using dual-tree ', num2str(sum(~angle_map_dt{lev}(:)))]); 
    display(['Number of mis-aligned points using ILP ', num2str(sum(~angle_map_ilp{lev}(:)))]);  
    display(['Number of mis-aligned points using ICP ', num2str(sum(~angle_map_icp{lev}(:)))]); 
end
%%
%Are the locations at which ICP aren't rotation invariant a subset of
%locations that have angles outside the maximal sub-bands...
all(abs(max_icp_diff_90{4}(:)) > 1e-3 == (abs(max_icp_diff_90{4}(:)) > 1e-3 & ~angle_map_icp{4}(:)))
all(abs(max_icp_diff_90{3}(:)) > 1e-3 == (abs(max_icp_diff_90{3}(:)) > 1e-3 & ~angle_map_icp{3}(:)))
all(abs(max_icp_diff_90{2}(:)) > 1e-3 == (abs(max_icp_diff_90{2}(:)) > 1e-3 & ~angle_map_icp{2}(:)))
all(abs(max_icp_diff_90{1}(:)) > 1e-3 == (abs(max_icp_diff_90{1}(:)) > 1e-3 & ~angle_map_icp{1}(:)))
%
all(abs(max_icp_diff_90{4}(:)) > 1e-3 == (abs(max_icp_diff_90{4}(:)) > 1e-3 & ~angle_map_dt{4}(:)))
all(abs(max_icp_diff_90{3}(:)) > 1e-3 == (abs(max_icp_diff_90{3}(:)) > 1e-3 & ~angle_map_dt{3}(:)))
all(abs(max_icp_diff_90{2}(:)) > 1e-3 == (abs(max_icp_diff_90{2}(:)) > 1e-3 & ~angle_map_dt{2}(:)))
all(abs(max_icp_diff_90{1}(:)) > 1e-3 == (abs(max_icp_diff_90{1}(:)) > 1e-3 & ~angle_map_dt{1}(:)))
% ... yes!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
for lev = 1:4

    figure('name', 'Dual-tree'); imagesc(max_icp_diff_90{lev} > 1e-3); axis image;
end
%%
%Does the problem only occur for small coefficients
for lev = 1:4

    [m_dt] = max(dt_lenna{lev}, [], 3);
    [m_ilp] = max(ilp_0{lev}, [], 3);
    [m_icp] = max(icp_0{lev}, [], 3);
 
    %Do difference occur in small or large magnitude coeffs
    display([mean(abs(m_dt(max_icp_diff_90{lev}(:) > 1e-3))), mean(abs(m_dt(max_icp_diff_90{4}(:) <= 1e-3)))]);
    display([mean(abs(m_icp(max_icp_diff_90{lev}(:) > 1e-3))), mean(abs(m_icp(max_icp_diff_90{4}(:) <= 1e-3)))]);
    display([mean(abs(m_ilp(max_icp_diff_90{lev}(:) > 1e-3))), mean(abs(m_ilp(max_icp_diff_90{4}(:) <= 1e-3)))]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%