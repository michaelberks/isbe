%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Looking again at the effects of different spline warps
%%
% Setting up source and target points

shapes_unaligned = get_shapes_from_masses(u_files1, 500);
shapes_aligned = align_shapes(shapes_unaligned, 'shiftOrigin', 1);
scaled_mean_shape = mass_model.mean_shape * mass_model.scale_factor;
%Define source points for TPS - as row vectors
s_x = scaled_mean_shape(1:500) + mass_model.mean_centre(1);
s_y = scaled_mean_shape(501:end) + mass_model.mean_centre(2);

%Define points to be interpolated by TPS - as row vectors
i_x = mass_model.mean_shape_pl(:,1)';
i_y = mass_model.mean_shape_pl(:,2)';

%Define displacement to target points
z_x = shapes_unaligned(1, 1:500);
z_y = shapes_unaligned(1, 501:end);

joint_plots = figure; hold on; axis equal;
%%
% Thin-plate
tps_L_inv = tps_weights(s_x, s_y, 'biharmTPS');

%Compute displacement of interpolated points
f_x = tps_warp(tps_L_inv, s_x, s_y, z_x, i_x, i_y, 'biharmTPS');
f_y = tps_warp(tps_L_inv, s_x, s_y, z_y, i_x, i_y, 'biharmTPS');

TPS = [f_x' f_y'];
figure; plot(TPS(:,1), TPS(:,2), 'r.'); axis equal;
figure(joint_plots); plot(TPS(:,1), TPS(:,2), 'r.');
%
% Clamped plate
tps_L_inv = tps_weights(s_x, s_y, 'biharmCPS');

%Compute displacement of interpolated points
f_x = tps_warp(tps_L_inv, s_x, s_y, z_x, i_x, i_y, 'biharmCPS');
f_y = tps_warp(tps_L_inv, s_x, s_y, z_y, i_x, i_y, 'biharmCPS');

CPS = [f_x' f_y'];
figure; plot(CPS(:,1), CPS(:,2), 'r.'); axis equal;
figure(joint_plots); plot(CPS(:,1), CPS(:,2), 'gx');
%%
% Now try geom_align
trans = geom_alignpoints_2d([z_x; z_y], [s_x; s_y], ones(1,500),...
    'transform', 'spline');
[TPS_isbe] = geom_transformpoints([i_x; i_y], trans);
TPS_isbe = TPS_isbe';

figure; plot(TPS_isbe(:,1), TPS_isbe(:,2), 'bo'); axis equal;
figure(joint_plots); plot(TPS_isbe(:,1), TPS_isbe(:,2), 'bo');
%%
trans = geom_alignpoints_2d([z_x; z_y], [s_x; s_y], ones(1,500),...
    'transform', 'spline', 'green', 'biharmCPS');
[CPS_isbe] = geom_transformpoints([i_x; i_y], trans);
CPS_isbe = CPS_isbe';

figure; plot(CPS_isbe(:,1), CPS_isbe(:,2), 'bo'); axis equal;
figure(joint_plots); plot(CPS_isbe(:,1), CPS_isbe(:,2), 'ys');

%%
% Now try geom_align
trans = geom_alignpoints_2d([z_x; z_y], [s_x; s_y], [],...
    'transform', 'spline');
[TPS_isbe2] = geom_transformpoints([i_x; i_y], trans);
TPS_isbe2 = TPS_isbe2';

figure; plot(TPS_isbe2(:,1), TPS_isbe2(:,2), 'm^'); axis equal;
figure(joint_plots); plot(TPS_isbe2(:,1), TPS_isbe2(:,2), 'm^');
%%
trans = geom_alignpoints_2d([z_x; z_y], [s_x; s_y], [],...
    'transform', 'spline', 'green', 'biharmCPS');
[CPS_isbe2] = geom_transformpoints([i_x; i_y], trans);
CPS_isbe2 = CPS_isbe2';

figure; plot(CPS_isbe2(:,1), CPS_isbe2(:,2), 'cv'); axis equal;
figure(joint_plots); plot(CPS_isbe2(:,1), CPS_isbe2(:,2), 'cv');