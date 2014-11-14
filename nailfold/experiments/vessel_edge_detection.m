load('C:\isbe\nailfold\data\pilot_study\apexes\normal\apex2872_vessel.mat')

g = gaussian_filters_1d(2);
g = g / sum(g);
v_pts_sm = [conv(v_pts(:,1), g, 'valid') conv(v_pts(:,2), g, 'valid')];


idx = 1:2:size(v_pts_sm,1);

[normal_xy, normal_ppx, normal_ppy, dists] = compute_spline_normals(v_pts_sm(idx,:));

figure; imgray(vessel_patch);
plot(v_pts_sm(:,1), v_pts_sm(:,2), '--');
plot(...
    [v_pts_sm(idx,1)+10*normal_xy(:,1) v_pts_sm(idx,1)-10*normal_xy(:,1)]',...
    [v_pts_sm(idx,2)+10*normal_xy(:,2) v_pts_sm(idx,2)-10*normal_xy(:,2)]', 'r');
plot(...
    v_pts_sm(idx,1)-10*normal_xy(:,1),...
    v_pts_sm(idx,2)-10*normal_xy(:,2), 'gx');
plot(...
    v_pts_sm(idx,1)+10*normal_xy(:,1),...
    v_pts_sm(idx,2)+10*normal_xy(:,2), 'yx');

%%
figure;
for i_n = 1:16
    norm_profile = improfile(apex_patch,...
        vessel_top(idx(i_n),1)+normal_xy(i_n,1)*[-30 30],...
        vessel_top(idx(i_n),2)+normal_xy(i_n,2)*[-30 30], 61);
        subplot(4,4,i_n); plot(11:51, norm_profile(11:51));
end
%%
sigma = 2;    
figure;
for i_n = 1:16
    norm_profile = improfile(apex_patch,...
        vessel_top(idx(i_n),1)+normal_xy(i_n,1)*[-30 30],...
        vessel_top(idx(i_n),2)+normal_xy(i_n,2)*[-30 30], 61); 
        
    [~, dg] = gaussian_filters_1d(sigma);
    normal_g1 = conv(norm_profile, dg, 'valid');
    
    norm_profile = (norm_profile - mean(norm_profile)) / (max(norm_profile) - min(norm_profile));
    %normal_g1 = normal_g1 / (max(normal_g1) - min(normal_g1));

    xo = (length(dg) -1)/2;
    x = (xo+1):(length(norm_profile)-xo);

    subplot(4,4,i_n); hold all;
    plot(1:61, norm_profile);
    plot(x, normal_g1);
    plot([31 31], [-1 1], 'k--');
    axis([1 61 -1 1]);

end
%%
sigma1 = 2;
sigma2 = 4;
x = 26:76;
idx_inner = 26:50;
idx_outer = 52:76;
figure; imgray(vessel_patch);
%plot(vessel_top(:,1), vessel_top(:,2));

num_pts = size(normal_xy,1);
outer_edge_xy = zeros(num_pts,2);
inner_edge_xy = zeros(num_pts,2);

for i_n = 1:num_pts
    [norm_xp norm_yp norm_profile] = improfile(vessel_patch,...
        v_pts_sm(idx(i_n),1)+normal_xy(i_n,1)*[-50 50],...
        v_pts_sm(idx(i_n),2)+normal_xy(i_n,2)*[-50 50], 101, 'bilinear');
    
    [~, dg] = gaussian_filters_1d(sigma1);
    normal_g1 = conv(norm_profile, dg, 'same');

    max_inner = find(diff(normal_g1(idx_inner)) <= 0, 1, 'last')+1;
    max_outer = find(diff(normal_g1(idx_outer)) <= 0, 1, 'first');
       
    x_edge_inner = norm_xp(idx_inner(max_inner));
    y_edge_inner = norm_yp(idx_inner(max_inner));
    x_edge_outer = norm_xp(idx_outer(max_outer));
    y_edge_outer = norm_yp(idx_outer(max_outer));
    
    plot(x_edge_inner, y_edge_inner, 'rx');
    plot(x_edge_outer, y_edge_outer, 'gx');
    
    outer_edge_xy(i_n,:) = [x_edge_outer y_edge_outer];
    inner_edge_xy(i_n,:) = [x_edge_inner y_edge_inner];

end
%%
edge_normal_xy = compute_spline_normals(outer_edge_xy);

match_edge_xy = size(inner_edge_xy);

for i_n = 1:num_pts
    ox = outer_edge_xy(i_n,1);
    oy = outer_edge_xy(i_n,2);
    
    inner_vecs = [inner_edge_xy(:,1)-ox inner_edge_xy(:,2)-oy];
    inner_dists = sum(inner_vecs.^2,2);
    inner_vecs = bsxfun(@rdivide, inner_vecs, sqrt(inner_dists));
    inner_dots = sum(inner_vecs .* -edge_normal_xy, 2);
    
    valid_idx = inner_dots > 0.9;
    [~, nearest_idx] = min(inner_dists(valid_idx));
    
    valid_edge_xy = inner_edge_xy(valid_idx,:);
    match_edge_xy(i_n,:) = valid_edge_xy(nearest_idx,:);
    
    display(inner_dots(nearest_idx));
end
figure; imgray(vessel_patch);
plot(...
    [outer_edge_xy(:,1) match_edge_xy(:,1)]',...
    [outer_edge_xy(:,2) match_edge_xy(:,2)]');
%%
med_inner_edge = medfilt1(inner_edge_xy, 5, [], 1);
%%
figure; axis equal ij; hold all;
plot(outer_edge_xy(:,1), outer_edge_xy(:,2), 'k');
plot(inner_edge_xy(:,1), inner_edge_xy(:,2));
plot(inner_edge_xy(:,1), inner_edge_xy(:,2), 'rx');
plot(med_inner_edge(:,1), med_inner_edge(:,2), 'c--');
plot(med_inner_edge(:,1), med_inner_edge(:,2), 'g.');

dist_mat = squareform(pdist(inner_edge_xy));
%%
med_inner_edge = medfilt1(match_edge_xy, 5, [], 1);
%
figure; axis equal ij; hold all;
plot(match_edge_xy(:,1), match_edge_xy(:,2));
plot(match_edge_xy(:,1), match_edge_xy(:,2), 'rx');
plot(med_inner_edge(:,1), med_inner_edge(:,2), 'c--');
plot(med_inner_edge(:,1), med_inner_edge(:,2), 'g.');
%%
prof_width = 20;
outer_rectangle_x = repmat(0:prof_width, num_pts, 1);
outer_rectangle_y = repmat((1:num_pts)', 1, prof_width+1);

matched_inner_xy = outer_edge_xy - prof_width*edge_normal_xy;


%Define source points for TPS - as row vectors
s_x = [outer_rectangle_x(:,1)' outer_rectangle_x(1:5:end,end)'];
s_y = [outer_rectangle_y(:,1)' outer_rectangle_y(1:5:end,end)'];

%Define points to be interpolated by TPS - as row vectors
i_x = outer_rectangle_x(:)';
i_y = outer_rectangle_y(:)';
   
%Define displacement to target points
z_x = [outer_edge_xy(:,1)' matched_inner_xy(1:5:end,1)'];
z_y = [outer_edge_xy(:,2)' matched_inner_xy(1:5:end,2)'];     

T = geom_alignpoints_2d([z_x; z_y], [s_x; s_y], [], 'transform', 'spline');
[pts] = geom_transformpoints([i_x; i_y], T);
%
figure; axis equal ij; hold all;
plot(i_x, i_y, 'r.');
plot(s_x, s_y, 'b.');

figure; axis equal ij; hold all;
plot(pts(1,:), pts(2,:), 'r.');
plot(z_x, z_y, 'b.');

%
vessel_pts_x = reshape(pts(1,:), [num_pts prof_width+1]);
vessel_pts_y = reshape(pts(2,:), [num_pts prof_width+1]);

figure; axis equal ij; hold all;
for i_p = 1:num_pts
    plot(vessel_pts_x(i_p,:), vessel_pts_y(i_p,:));
end
plot(outer_edge_xy(:,1), outer_edge_xy(:,2), 'k');
%%
figure; axis equal ij; hold all;
plot([outer_edge_xy(:,1) inner_edge_xy(:,1)]', [outer_edge_xy(:,2) inner_edge_xy(:,2)]');
%%
oe = outer_edge_xy;
ie = inner_edge_xy;
normals_cross = false(num_pts, num_pts);
for i_p = 1: num_pts - 10
    
    vec1 = [oe(i_p,:); ie(i_p,:)];
    for j_p = i_p + (1:10)
        %figure; axis equal ij; hold all;
        %plot(vec1(:,1), vec1(:,2));
        %plot(vec2(:,1), vec2(:,2));
        
        vec2 = [oe(j_p,:); ie(j_p,:)];
        
        normals_cross(i_p, j_p) = line_cross(vec1, vec2);
        normals_cross(j_p, i_p) = normals_cross(i_p, j_p);
    end
end
figure; imgray(normals_cross);
%%
keep_idx = (1:num_pts)';
while any(normals_cross(:))
    [~, remove_i] = max(sum(normals_cross));
    display(remove_i);
    
    normals_cross(remove_i,:) = [];
    normals_cross(:,remove_i) = [];
    keep_idx(remove_i) = [];
end
%%
figure; axis equal ij; hold all;
plot(...
    [outer_edge_xy(keep_idx,1) inner_edge_xy(keep_idx,1)]',...
    [outer_edge_xy(keep_idx,2) inner_edge_xy(keep_idx,2)]');
%%
[outer_edge_xy, inner_edge_xy, vessel_mask] = detect_vessel_edge(vessel_patch, v_pts);
%%
v_files = dir('C:\isbe\nailfold\data\pilot_study\apexes\normal\*vessel.mat');
for i_v = 151:200
    load(['C:\isbe\nailfold\data\pilot_study\apexes\normal\' v_files(i_v).name]);
    
    if strcmpi(apex_properties.vessel_shape, 'Normal')
        detect_vessel_edge(vessel_patch, v_pts);
    end
end
%%
s1 = load('C:\isbe\nailfold\data\pilot_study\apexes\normal\apex2872.mat');
s2 = load('C:\isbe\nailfold\data\pilot_study\apexes\normal\apex2872_vessel.mat');

col_offset = s1.apex_properties.sc - s2.vessel_properties.sc;
row_offset = s1.apex_properties.sr - s2.vessel_properties.sr;

vessel_patch = double(s2.vessel_patch);

g = gaussian_filters_1d(16, 48);
g = g / sum(g);
im_edges = conv2(g', g, ones(size(vessel_patch)), 'same');
vessel_patch_smoothed = conv2(g', g, vessel_patch, 'same') ./ im_edges;
vessel_patch_equalised = vessel_patch - vessel_patch_smoothed;

v_pts = s2.v_pts;
apex_xy = [s1.apex_xy(:,1)+col_offset, s1.apex_xy(:,2)+row_offset];

diffs = (v_pts(:,1) - mean(apex_xy(:,1))).^2 + (v_pts(:,2) - mean(apex_xy(:,2))).^2;
[~, nearest_pt] = min(diffs);

[~,~,nearest_pt2] = line_cross(v_pts, apex_xy);

figure; imgray(vessel_patch_equalised);
plot(apex_xy(:,1), apex_xy(:,2), 'r-x');
plot(v_pts(:,1), v_pts(:,2), 'g');
plot(v_pts(:,1), v_pts(:,2), 'g.');
plot(v_pts(nearest_pt,1), v_pts(nearest_pt,2), 'gx');
plot(v_pts(nearest_pt2,1), v_pts(nearest_pt2,2), 'cx');
%%                
detect_vessel_edge_snake(vessel_patch_equalised, v_pts, apex_xy);
%%
for i_size = [0.5 1 2 4]
    
    vessel_patch_l = imresize(vessel_patch, i_size, 'bicubic');
    apex_vec = diff(wide_apex_xy);
    theta = atan(apex_vec(2)/apex_vec(1));
    cc = cos(theta).^2;
    ss = sin(theta).^2;
    s2 = sin(2*theta);

    sigma = [0.01 1 2 4 8 16 32];

    [g2d_responses] = compute_gaussian_2nd_derivatives(vessel_patch_l, sigma);


    figure; hold all;
    for i_scale = 1:length(sigma)

        Ixy = improfile(g2d_responses(:,:,i_scale,1), i_size*wide_apex_xy(:,1), i_size*wide_apex_xy(:,2), 'bicubic');
        Ixx = improfile(g2d_responses(:,:,i_scale,2), i_size*wide_apex_xy(:,1), i_size*wide_apex_xy(:,2), 'bicubic');
        Iyy = improfile(g2d_responses(:,:,i_scale,3), i_size*wide_apex_xy(:,1), i_size*wide_apex_xy(:,2), 'bicubic');


        wo_theta = Ixx.*cc + Iyy.*ss + Ixy.*s2;
        plot(wo_theta);
        %figure; imgray(smooth_patch);
    end
end
%%
for sigma = [1 2 4 8 16]
    %smooth_patch = imfilter(vessel_patch, fspecial('gaussian', 6*sigma, sigma), 'replicate');
    g = gaussian_filters_1d(sigma);
    g = g / sum(g);
    g2 = g'*g;
    figure; imgray(conv2(g',g,ones(128))); colorbar;
end
%%
apex_width = sqrt(sum(diff(apex_xy).^2));
for i_scale = [1 2 4 8]
    [g2d_responses] = compute_gaussian_2nd_derivatives(vessel_patch, apex_width/i_scale);
    [g2d_mag, g2d_ori] = gaussian_2nd_derivative_line(g2d_responses);
    g2d_combined = g2d_mag.*complex(cos(2*g2d_ori), sin(2*g2d_ori));
    figure; imgray(complex2rgb(g2d_combined));
    nms_vessels = mb_non_maximal_supp(g2d_mag,g2d_ori);
    figure; imgray(nms_vessels > 0);
end
    %figure; mesh(g2d_mag);
%%
f1 = figure;
display_orientation(f1, g2d_mag, complex(cos(g2d_ori), sin(g2d_ori)), 2);
%%
[g2d_responses] = compute_gaussian_2nd_derivatives(vessel_patch, apex_width/2);
[g2d_mag, g2d_ori] = gaussian_2nd_derivative_line(g2d_responses);
%%
x_grid = 1:2:129;
y_grid = (1:2:291)';
for mu = (1)/100
    [u,v] = GVF(g2d_mag, mu, 1000);
    figure; axis equal ij; hold on; quiver(x_grid, y_grid, u(1:2:end,1:2:end), v(1:2:end,1:2:end));
    plot(v_pts(:,1), v_pts(:,2), 'r');
end
%%
x_grid = 1:2:129;
y_grid = (1:2:291)';
x = v_pts(:,1);
y = v_pts(:,2);
x = [x(1); x(1); x; x(end); x(end)];
y = [y(1); y(1); y; y(end); y(end)];
alpha = 0;
beta = 0;
gamma = 1;
kappa = 0.6;

[px,py] = GVF(g2d_mag, 0.2, 100);
for i=1:20,
    [x,y] = snakedeform(x, y, alpha, beta, gamma, kappa, px, py, 20);
    %[x,y] = snakeinterp1(x,y,2,0.5);
    figure; axis equal ij; hold on;
    quiver(x_grid, y_grid, u(1:2:end,1:2:end), v(1:2:end,1:2:end));
    plot(v_pts(:,1), v_pts(:,2), 'r');
    plot(x, y, 'g');
    title(['Deformation in progress,  iter = ' num2str(i*5)])
end
x = x(3:end-2);
y = y(3:end-2);

figure; imgray(vessel_patch_equalised);
plot(x, y, 'g');
%%
v_files = dir('C:\isbe\nailfold\data\pilot_study\apexes\normal\*vessel.mat');

for i_v = 700%:160
    s1 = load(['C:\isbe\nailfold\data\pilot_study\apexes\normal\' v_files(i_v).name(1:8) '.mat']);
    s2 = load(['C:\isbe\nailfold\data\pilot_study\apexes\normal\' v_files(i_v).name]);

    col_offset = s1.apex_properties.sc - s2.vessel_properties.sc;
    row_offset = s1.apex_properties.sr - s2.vessel_properties.sr;

    vessel_patch = double(s2.vessel_patch);

    g = gaussian_filters_1d(16, 48);
    g = g / sum(g);
    im_edges = conv2(g', g, ones(size(vessel_patch)), 'same');
    vessel_patch_smoothed = conv2(g', g, vessel_patch, 'same') ./ im_edges;
    vessel_patch_equalised = vessel_patch - vessel_patch_smoothed;

    v_pts = s2.v_pts;
    apex_xy = [s1.apex_xy(:,1)+col_offset, s1.apex_xy(:,2)+row_offset];
    try          
        detect_vessel_edge(vessel_patch_equalised, v_pts, apex_xy, 'plot', 1, 'initial_spacing', 4);
    catch last_err
        display(['Error in vessel ' num2str(i_v)]);
        display(last_err.message);
    end
end
%%

v_files = dir('C:\isbe\nailfold\data\rsa_study\apexes\normal\*vessel.mat');
for i_v = 81:100
    s = load(['C:\isbe\nailfold\data\rsa_study\apexes\normal\' v_files(i_v).name]);
    
    vessel_patch = double(s.vessel_patch);

%     g = gaussian_filters_1d(16, 48);
%     g = g / sum(g);
%     im_edges = conv2(g', g, ones(size(vessel_patch)), 'same');
%     vessel_patch_smoothed = conv2(g', g, vessel_patch, 'same') ./ im_edges;
%     vessel_patch_equalised = vessel_patch - vessel_patch_smoothed;
    
    figure; imgray(vessel_patch);
    plot(s.v_pts(:,1), s.v_pts(:,2));
end
