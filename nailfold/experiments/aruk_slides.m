%%
load([nailfoldroot 'data\apex_detection\apex_templates_aligned.mat'], 'apex_template*');
load([nailfoldroot 'data\apex_detection\aligned\thresh.mat'], 'thresh2');
load([nailfoldroot 'data\apex_detection\mean_shape.mat'], 'mean_shape');

pilot_dir = 'C:\isbe\nailfold\data\pilot_study\';
%%
image_paths = {...
    ['C:\isbe\nailfold\images\camera_demos\2012_08_06\mike_mosaic\'...
        'Visit1\Lhand\Digit4\X300\FullResMosaic\STUDENT2V1LD4X3LrgMosaic.bmp'];...
    ['C:\isbe\nailfold\images\camera_demos\2012_08_06\mike_mosaic\'...
        'Visit1\Lhand\Digit5\X300\FullResMosaic\STUDENT2V1LD5X3LrgMosaic.bmp'];...
    [pilot_dir 'images\14735c.png'];...
    [pilot_dir 'images\40013c.png'];...
    [pilot_dir 'images\59214c.png'];...
    'C:\isbe\nailfold\data\rsa_study\images\10598c.png'};

aam_paths = {...
	'C:\isbe\nailfold\data\d4\candidates\template_matching\aligned\';...
    'C:\isbe\nailfold\data\d5\candidates\template_matching\aligned\';...
    'C:\isbe\nailfold\data\p14735\candidates\template_matching\aligned\';...
    'C:\isbe\nailfold\data\p40013\candidates\template_matching\aligned\';...
    'C:\isbe\nailfold\data\p59214\candidates\template_matching\aligned\';...
    'C:\isbe\nailfold\data\p10598\candidates\template_matching\aligned\'};
    

%%
for i_im = 6
    detect_capillaries(image_paths{i_im}, ...
        'g2d_template', apex_template2_a,...
        'template_thresh', thresh2,...
        'max_num_candidates', 50,...
        'thetas', -15:3:15,...
        'scales', 0.8:0.1:1.5,...
        'mean_shape', mean_shape,...
        'aam_dir', aam_paths{i_im},...
        'aam_exe', 'ncm_sandpit_mb',...
        'aam_path', 'C:\isbe\nailfold\playground\model_experiments\models\orig\2\vessel_apex_orig.smd',...
        'g2d_vals', [],...
        'g2d_thresh', 0.9,...
        'do_template_matching', 1,...
        'do_aam', 1, ...
        'do_final_vessel', 1, ...
        'delete_candidate_patches', 1,...
        'nailfold_mask', []);
end
%%
for i_im = 1:2
    nailfold = imread(image_paths{i_im});
    mask = make_nailfold_mosaic_mask(nailfold, 240, 20);

    nailfold_r = double(nailfold);
    clear nailfold;
    min_g = min(nailfold_r(mask));
    max_g = max(nailfold_r(mask));

    nailfold_r = (nailfold_r - min_g) / (max_g - min_g);
    nailfold_r(nailfold_r < 0) = 0;
    nailfold_r(nailfold_r > 1) = 1;

    nailfold_g = nailfold_r;
    nailfold_b = nailfold_r;

    capillary_colors = hsv(50);

    for i_ap = 1:50
        load([aam_paths{i_im} 'apex' zerostr(i_ap, 4) '_aam.mat'], 'apex_candidate');
        if isfield(apex_candidate, 'is_distal') && apex_candidate.is_distal
            apex_x = apex_candidate.final_vessel(:,1) + apex_candidate.sc;
            apex_y = apex_candidate.final_vessel(:,2) + apex_candidate.sr;
            apex_idx = sub2ind(size(mask), apex_y, apex_x);

            nailfold_r(apex_idx) = capillary_colors(i_ap, 1);
            nailfold_g(apex_idx) = capillary_colors(i_ap, 2);
            nailfold_b(apex_idx) = capillary_colors(i_ap, 3);
        end
    end

    nailfold_im = cat(3, nailfold_r, nailfold_g, nailfold_b);

    figure; imgray(nailfold_im);
    imwrite(nailfold_im, ['im' num2str(i_im) '.png']);
end

%%
i_im = 1;
i_ap = 1;
nailfold = imread(image_paths{i_im});
load([aam_paths{i_im} 'apex' zerostr(i_ap, 4) '_aam.mat'], 'apex_candidate');

nailfold_patch = nailfold(apex_candidate.sr:apex_candidate.er, apex_candidate.sc:apex_candidate.ec);
vessel_idx = sub2ind(size(nailfold_patch), apex_candidate.final_vessel(:,2), apex_candidate.final_vessel(:,1));
vessel_patch = false(size(nailfold_patch));
vessel_patch(vessel_idx) = 1;

perim_px = bwperim(vessel_patch);
[perim_y perim_x] = find(perim_px);
%%

x = apex_candidate.fitted_vessel_xy(:,1);
y = apex_candidate.fitted_vessel_xy(:,2);
d = cumsum([0; sum(diff([x y]).^2,2)]);

ppx = spline(d, x);
ppy = spline(d, y);

ppx2 = ppx;
ppx2.order = 3;
ppx2.coefs = ppx.coefs(:,1:3)*diag([3 2 1]);

ppy2 = ppy;
ppy2.order = 3;
ppy2.coefs = ppy.coefs(:,1:3)*diag([3 2 1]);

yg = ppval(ppy2, d);
xg = ppval(ppx2, d);
norm_x = -yg ./ sqrt(xg.^2 + yg.^2);
norm_y = xg ./ sqrt(xg.^2 + yg.^2);
%
figure; imgray(nailfold_patch);

for i_p = [1 16 31]
    perim_xd = perim_x - x(i_p); 
    perim_yd = perim_y - y(i_p);
    perim_l = sqrt(perim_xd.^2 + perim_yd.^2);
    perim_xd = perim_xd ./ perim_l;
    perim_yd = perim_yd ./ perim_l;
    
    dot_prods = perim_xd*norm_x(i_p) + perim_yd*norm_y(i_p);
    inner_idx = dot_prods > 0.5;
    outer_idx = dot_prods < -0.5;
    
    inner_perim_x = perim_x(inner_idx);
    inner_perim_y = perim_y(inner_idx);
    outer_perim_x = perim_x(outer_idx);
    outer_perim_y = perim_y(outer_idx);
    
    [~,inner_min] = min(perim_l(inner_idx));
    [~,outer_min] = min(perim_l(outer_idx));
    
    edge_x1 = inner_perim_x(inner_min);
    edge_y1 = inner_perim_y(inner_min);
    edge_x2 = outer_perim_x(outer_min);
    edge_y2 = outer_perim_y(outer_min);
    
    plot([edge_x1 edge_x2], [edge_y1 edge_y2], 'r-x', 'linewidth', 3);
end
axis off;
saveas(gcf, 'apex_width.png');
%%
figure; imgray(nailfold_patch);
plot(x+1, y+1, 'gx', 'markersize', 10);
axis off;
saveas(gcf, 'apex_aam.png');

write_im_from_colormap(nailfold_patch, 'apex.png');
% plot([x + 10*norm_x, x - 10*norm_x]', [y + 10*norm_y y - 10*norm_y]');
%% 
widths = [];
for i_im = 1
    for i_ap = 1:50
        load([aam_paths{i_im} 'apex' zerostr(i_ap, 4) '_aam.mat'], 'apex_candidate');
        if isfield(apex_candidate, 'is_distal') && apex_candidate.is_distal
            widths(end+1) = apex_candidate.mean_apex_widths;
        end
    end
end    
%% 
apex_centres = zeros(0,2);
for i_im = 1
    for i_ap = 1:50
        load([aam_paths{i_im} 'apex' zerostr(i_ap, 4) '_aam.mat'], 'apex_candidate');
        if isfield(apex_candidate, 'is_distal') && apex_candidate.is_distal
            x = apex_candidate.fitted_vessel_xy(16,1) + apex_candidate.sc;
            y = apex_candidate.fitted_vessel_xy(16,2) + apex_candidate.sr;
            apex_centres(end+1,:) = [x y];
        end
    end
end    
apex_centres = sortrows(apex_centres);
figure; imgray(nailfold);
plot(apex_centres(:,1), apex_centres(:,2));

mean(diff(apex_centres))*1.183
dists = sqrt(sum(diff(apex_centres).^2,2));
    

%%
figure; imgray(perim_px);
plot([x + 10*norm_x, x - 10*norm_x]', [y + 10*norm_y y - 10*norm_y]');
plot(x, y, 'r', 'linewidth', 3);
%%
nailfold = imread(image_path5);
mask = make_nailfold_mosaic_mask(nailfold, 240, 20);

figure; imgray(nailfold); caxis([min(nailfold(mask)) max(nailfold(mask))]);

%
for i_ap = 1:50
    load([aam_dir5 'apex' zerostr(i_ap, 4) '_aam.mat'], 'apex_candidate');
    apex_x = apex_candidate.fitted_vessel_xy(:,1) + apex_candidate.sc;
    apex_y = apex_candidate.fitted_vessel_xy(:,2) + apex_candidate.sr;
    
    plot(apex_x, apex_y);
    text(apex_x(end), apex_y(end), num2str(i_ap), 'color', 'r');
    
    if i_ap == 1
        current_apexes = [apex_x(16) apex_y(16)];
    else
        tan_vectors = abs(atan((current_apexes(:,2)-apex_y(16)) ./ (current_apexes(:,1)-apex_x(16))));
        if any(tan_vectors > pi/4)
            plot(apex_x(16), apex_y(16), 'rx', 'markersize', 10);
        else
            current_apexes = [current_apexes; apex_x(16) apex_y(16)]; %#ok
            plot(apex_x(16), apex_y(16), 'gx', 'markersize', 10);
        end
    end   
end
%%
for i_ap = 1:32
    load([aam_dir 'apex' zerostr(i_ap, 4) '_aam.mat'], 'apex_candidate');
       
    nailfold_patch = im4(...
        apex_candidate.sr:apex_candidate.er,...
        apex_candidate.sc:apex_candidate.ec);
    
    figure; imgray(nailfold_patch);   
    plot(apex_candidate.fitted_vessel_xy(:,1), apex_candidate.fitted_vessel_xy(:,2));
    title(['Model score: ' num2str(apex_candidate.model_score)]);
    
    if ~isfield(apex_candidate, 'final_vessel')
        display([num2str(i_ap) ' missing']);
        continue;
    end
    
    plot(apex_candidate.final_vessel(:,1), apex_candidate.final_vessel(:,2), 'rx');
    
end
    
    
%%
pilot_dir = 'C:\isbe\nailfold\data\pilot_study\';
load([pilot_dir '\models\apex_templates.mat'], 'mean_apex_*');

im4 = imread('C:\isbe\nailfold\images\camera_demos\2012_08_06\mike_mosaic\Visit1\Lhand\Digit4\X300\FullResMosaic\STUDENT2V1LD4X3LrgMosaic.bmp');
mask = make_nailfold_mosaic_mask(im4, 240, 20);
%%
templates = cell(2,2);
templates{1,1} = 'g1d';
templates{2,1} = 'g2d';

%Original
templates{1,2} = mean_apex_g1(:,:,1,1);
templates{2,2} = mean_apex_g2(:,:,1,1);
            
[maxima_pos, maxima_vals] = template_match_apexes(im4, templates, 'plot', 0);

%%
figure; imgray(im4);
plot(maxima_pos(1:40,1), maxima_pos(1:40,2), 'rx');
%%
[g2d_map, ori_map] = gaussian_2nd_derivative_line(im4, 4);

nms_vessels = mb_non_maximal_supp(g2d_map, ori_map);
nms_vessels_bw = bwareaopen( nms_vessels > 0, 20);

figure; imgray(nms_vessels);

x = repmat(-10:10, 21, 1);
y = x';
circ_pts = (x.^2 + y.^2) < 10^2;
circ_x = x(circ_pts);
circ_y = y(circ_pts);
%%
connected_vessels = false(size(im4));
for i_v = 1:40
    connected_vessels = connected_vessels | bwselect(nms_vessels_bw, maxima_pos(i_v,1)+circ_x, maxima_pos(i_v,2)+circ_y);
end

figure; imgray(connected_vessels);
plot(maxima_pos(1:40,1), maxima_pos(1:40,2), 'rx');    
%%
vals = sort(g2d_map(mask));
prctiles = linspace(0,1,length(vals));

g2d_probs = interp1(vals, prctiles, g2d_map);
g2d_probs(~mask) = 0;
figure; imgray(g2d_probs);
%%
poss_vessels = g2d_probs > .9;
connected_vessels = false(size(im4));
for i_v = 1:40
     vessels_i = bwselect(poss_vessels, maxima_pos(i_v,1), maxima_pos(i_v,2));
     vessels_is = bwmorph(vessels_i, 'skel', inf);
     vessels_is = bwareaopen(vessels_is, 100);
     
     [sy sx] = find(vessels_is);
     vessels_i = bwselect(vessels_i, sx, sy);
     
     connected_vessels = connected_vessels | vessels_i;
end
figure; imgray(connected_vessels);
%%
all_paths = zeros(size(im4));
I_ori_D = 0.8*ones(size(im4));
for i_v = 1:40
    for i_step = 1:1e3
        path_map_i = prob_track(g2d_probs, ori_map, I_ori_D, maxima_pos(i_v,1), maxima_pos(i_v,2), 2);
        all_paths = all_paths + path_map_i;
    end
    figure; imgray(all_paths);
end
%%
%--------------------------------------------------------------------------
%% 3) Create regions to test AAM fitting after initial candidate point selection
%--------------------------------------------------------------------------
%Reload threshold for aligned method
load([nailfoldroot 'data\apex_detection\apex_templates_aligned.mat'], 'apex_template*');
load([nailfoldroot 'data\apex_detection\aligned\thresh.mat'], 'thresh2');
load([nailfoldroot 'data\apex_detection\mean_shape.mat'], 'mean_shape');

x = repmat(-24:24, 49, 1);
y = repmat(-24:24, 49, 1)';
xy = [x(:) y(:)];

circle_mask = ~(x.^2 + y.^2 > 24^2);
N = sum(circle_mask(:));
template = apex_template2_a(circle_mask);

T = sum(template) / N;
T2 = sum(template.^2) / N;
denom_T = sqrt( max(T2 - T^2,0) );

apex_num = 1;
%             
%load image
image_path = 'C:\isbe\nailfold\images\camera_demos\2012_08_06\mike_mosaic\Visit1\Lhand\Digit4\X300\FullResMosaic\STUDENT2V1LD4X3LrgMosaic.bmp';
nailfold = imread(image_path);   
nailfold = nailfold(:,:,1);
nailfold_mask = imerode(nailfold > 10 & nailfold < 250, strel('disk', 20));

%Compute the gaussian derivatives
[mag_2d] = gaussian_2nd_derivative_line(nailfold, 4);

%Apply template matching
C2 = mb_normxcorr2(apex_template2_a, mag_2d);

%Get local maxima using precomputed threshold
[maxima_pos, maxima_vals] = local_image_maxima(C2, 20, nailfold_mask, thresh2);

num_maxima = size(maxima_pos,1);
maxima_scales = zeros(num_maxima,1);
maxima_thetas = zeros(num_maxima,1);

mkdir('C:\isbe\nailfold\data\fred\candidates\template_matching\aligned\');
for ii = 1:num_maxima
    c_max = 0;
    for theta = -15:3:15
        rot = [cosd(theta) sind(theta); -sind(theta) cosd(theta)];

        for scale = 0.8:0.1:1.5

            xya = xy * rot * scale;
            xa = reshape(xya(:,1) + maxima_pos(ii,1), 49, 49);
            ya = reshape(xya(:,2) + maxima_pos(ii,2), 49, 49);

            gd2_patch = interp2(mag_2d, xa, ya, 'bilinear');
            gd2_vec = gd2_patch(circle_mask);

            A = sum(gd2_vec) / N;
            A2 = sum(gd2_vec.^2) / N;

            denom = denom_T * sqrt( max(A2 - A^2,0) );
            numerator = (template' * gd2_vec) / N - A*T;

            c = numerator / denom;

            if c > c_max
                maxima_scales(ii) = scale;
                maxima_thetas(ii) = theta;
                vessel_xy = mean_shape * rot * scale;
                c_max = c;
            end

            %figure; imgray(gd2_patch);
            %title(['\theta = ' num2str(theta) ', scale = ' num2str(scale) ', C = ' num2str(c)]);

        end
    end

    sr = max(1, floor(min(vessel_xy(:,2)) - 50 + maxima_pos(ii,2)));
    er = min(size(nailfold,1), ceil(max(vessel_xy(:,2)) + 50 + maxima_pos(ii,2)));
    sc = max(1, floor(min(vessel_xy(:,1)) - 50 + maxima_pos(ii,1)));
    ec = min(size(nailfold,2), floor(max(vessel_xy(:,1)) + 50 + maxima_pos(ii,1)));

    vessel_xy(:,1) = vessel_xy(:,1) - sc + maxima_pos(ii,1);
    vessel_xy(:,2) = vessel_xy(:,2) - sr + maxima_pos(ii,2);

    %Sample patch from image
    image_patch = nailfold(sr:er, sc:ec);

    apex_candidate.image_patch = image_patch;
    apex_candidate.vessel_xy = vessel_xy;
    apex_candidate.true_vessel_xy = true_vessel_xy;
    apex_candidate.scale = maxima_scales(ii);
    apex_candidate.theta = maxima_thetas(ii);
    apex_candidate.sr = sr;
    apex_candidate.sc = sc;
    apex_candidate.method = 'template_matching_aligned_g2d';
    apex_candidate.image_path = image_path;

    save(['C:\isbe\nailfold\data\fred\candidates\template_matching\aligned\apex' zerostr(apex_num, 4) '.mat'], 'apex_candidate');
    apex_num = apex_num + 1;

end
%%
aam_dir = 'C:\isbe\nailfold\data\fred\candidates\template_matching\aligned\';
num_pts = 31;

fid(1) = fopen([aam_dir 'all_candidates_test.smd'], 'wt');

ff = 1;
fprintf(fid(ff), '%s \n', '// Text file describing test data for vessel apex models');
fprintf(fid(ff), '\n');
fprintf(fid(ff), '%s \n', '// Directory containing images');
fprintf(fid(ff), '%s \n', ['image_dir: ' aam_dir 'images/orig/']);
fprintf(fid(ff), '%s \n', '// Directory containing points');
fprintf(fid(ff), '%s \n', 'points_dir: C:/isbe/nailfold/data/aam/candidates/template_matching/aligned/points/');
fprintf(fid(ff), '\n');
fprintf(fid(ff), '%s \n', '// Details of points : images');
fprintf(fid(ff), '\n');
fprintf(fid(ff), '%s \n', 'training_set:');
fprintf(fid(ff), '%s \n', '{');

mkdir([aam_dir '\images\orig\']);
mkdir([aam_dir '\points\']);
for ii = 1:120
    
    %Load in saved candidate details
    load([aam_dir 'apex' zerostr(ii, 4) '.mat'], 'apex_candidate');
    vessel_xy = apex_candidate.vessel_xy;
    image_patch = apex_candidate.image_patch;
    clear apex_candidate;
        
    %Write out image patch
    imwrite(image_patch, [aam_dir '\images\orig\candidate_apex' zerostr(ii,4) '.png']);
    
    %Write out a pts file we can read in to VXL
    fid1 = fopen([aam_dir 'points\candidate_apex' zerostr(ii,4) '.pts'], 'wt');
    fprintf(fid1, '%s \n', 'version: 1');
    fprintf(fid1, '%s %d \n', 'n_points:', num_pts);
    fprintf(fid1, '%s \n', '{'); 
    for jj = 1:num_pts
        fprintf(fid1,'%.2f %.2f \n', vessel_xy(jj,1), vessel_xy(jj,2));
    end
    fprintf(fid1, '%s \n', '}');
    fprintf(fid1, 'nailfold: %s \n', image_path);
    fprintf(fid1, '%s %d \n', 'start_row:', sr);
    fprintf(fid1, '%s %d \n', 'start_col: ', sc);
    fclose(fid1);
    
    %Write entry for this image/pts pair in model .smd files
    str = ['candidate_apex' zerostr(ii,4) '.pts : candidate_apex' zerostr(ii,4) '.png'];
    fprintf(fid(1), '%s \n', str);
end

fprintf(fid(ff), '%s \n', '}');
fclose(fid(ff));
%%
fid = fopen([aam_dir 'out_points\model_qualities.txt']);
q_txt = textscan(fid, '%s %f', 'delimiter', ':');
fclose(fid);
model_q = q_txt{2}; clear q_txt;
[sorted_model_qualities qidx] = sort(model_q, 'descend');

feature = 'orig';

num_rows = 3;
num_cols = 4;
ii = 1;

while ii <= 32
    figure;
    for row = 1:num_rows, 
        for col = 1:num_cols
            if ii > 32
                break;
            end
            jj = qidx(ii);
            
            load([aam_dir 'apex' zerostr(jj, 4) '.mat'], 'apex_candidate');
            f1 = fopen([aam_dir 'out_points\candidate_apex' zerostr(jj,4) '.pts']);
            textscan(f1, '%[^{]');
            fgetl(f1);
            vessel_str = textscan(f1,'%[^}]');
            fclose(f1);
            test_pts = str2num(vessel_str{1}{1});

            im = imread([aam_dir 'images\orig\candidate_apex' zerostr(jj,4) '.png']);

            axes('units', 'normalized', 'position', [(col-1)/num_cols (num_rows-row)/num_rows 1/num_cols 1/num_rows]);
            
            [r c] = size(im);
            roi = im(round(r/4):round(3*r/4), round(c/4):round(3*c/4));
            
            imagesc(im); axis image off; colormap(gray(256)); caxis([min(roi(:)) max(roi(:))]);
            hold on;

            %Adjust the test points to be evenly spread
            dists = cumsum([0; sum(diff(test_pts).^2,2)]);
            test_pts = interp1(dists, test_pts, linspace(0, dists(end), 31));
            
            %plot(apex_candidate.vessel_xy(:,1),apex_candidate.vessel_xy(:,2), 'rx'); 
            plot(test_pts(:,1), test_pts(:,2), 'gx');                            
            
            %get normal vector at centre point
            x_norm = mean(test_pts(16:17,2) - test_pts(15:16,2));
            y_norm = mean(test_pts(15:16,1) - test_pts(16:17,1));
            plot(test_pts(16,1)+10*[-x_norm x_norm], test_pts(16,2)+10*[-y_norm y_norm], 'b-x');
            plot(test_pts(16,1), test_pts(16,2), 'rx');
            
            text(10, 10, sprintf('QoF: %5.2f', sorted_model_qualities(ii)), 'color', 'r');
            ii = ii + 1;
            
        end
    end
    
end
%%
image_path = 'C:\isbe\nailfold\images\camera_demos\2012_08_06\mike_mosaic\Visit1\Lhand\Digit4\X300\FullResMosaic\STUDENT2V1LD4X3LrgMosaic.bmp';
nailfold = imread(image_path);   
nailfold = nailfold(:,:,1);

fid = fopen([aam_dir 'out_points\model_qualities.txt']);
q_txt = textscan(fid, '%s %f', 'delimiter', ':');
fclose(fid);
model_q = q_txt{2}; clear q_txt;
[sorted_model_qualities qidx] = sort(model_q, 'descend');

feature = 'orig';

num_rows = 3;
num_cols = 4;
ii = 1;

figure; imgray(nailfold);
for ii = 1:32
    jj = qidx(ii);

    load([aam_dir 'apex' zerostr(jj, 4) '.mat'], 'apex_candidate');
    f1 = fopen([aam_dir 'out_points\candidate_apex' zerostr(jj,4) '.pts']);
    textscan(f1, '%[^{]');
    fgetl(f1);
    vessel_str = textscan(f1,'%[^}]');
    fclose(f1);
    test_pts = str2num(vessel_str{1}{1});
    test_pts(:,1) = test_pts(:,1) + apex_candidate.sc;
    test_pts(:,2) = test_pts(:,2) + apex_candidate.sr;

    plot(test_pts(:,1), test_pts(:,2), 'gx');
    text(test_pts(end,1), test_pts(end,2), num2str(ii), 'color', 'r');
    
end