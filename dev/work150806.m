p11 = 3.348;
p01 = *(pixel-istep);
  const T& p21 = *(pixel+istep);
  const T& p10 = *(pixel-jstep);
  const T& p12 = *(pixel+jstep);
  const T& p00 = *(&p10-istep);
  const T& p20 = *(&p10+istep);
  const T& p02 = *(&p12-istep);
  const T& p22 = *(&p12+istep);
%%
patch = v(83:85,142:144);
p00 = patch(1,1);
p10 = patch(1,2);
p20 = patch(1,3);

p01 = patch(2,1);
p11 = patch(2,2);
p21 = patch(2,3);

p02 = patch(3,1);
p12 = patch(3,2);
p22 = patch(3,3);

Ix =(-p00-p01-p02 +p20+p21+p22)/6.0;
Iy =(-p00-p10-p20 +p02+p12+p22)/6.0;
Ixx = ((p00+p01+p02 +p20+p21+p22)-2.0*(p10+p11+p12))/3.0;
Ixy = (p00+p22 -p02-p20)/4.0;
Iyy = ((p00+p10+p20 +p02+p12+p22)-2.0*(p01+p11+p21))/3.0;

det = Ixx*Iyy - Ixy*Ixy;
dx = (Iy*Ixy - Ix*Iyy) / det;
dy = (Ix*Ixy - Iy*Ixx) / det;
Io = p11;

val = Io + (Ix + 0.5*Ixx*dx + Ixy*dy)*dx + (Iy + 0.5*Iyy*dy)*dy;
%%
x = linspace(1,3,100);
y = linspace(1,3,100);
q = interp2(patch, x, y', '*bicubic');
figure; mesh(q);
[m i] = max(q(:))
[r c] = ind2sub([100 100], i)
[x(c) y(r)] - 2;
%%
dxx = repmat(x, 100, 1) - 2;
dyy = repmat(y', 1, 100) - 2;
q2 = Io + (Ix + 0.5*Ixx*dxx + Ixy*dyy).*dxx + (Iy + 0.5*Iyy*dyy).*dyy;
figure; mesh(q2);
%%
im = imread('C:\isbe\nailfold\data\rsa_study\cxx2\images\enlargedapex0140_vessel.png');
detect_capillaries(im, 1:4, ...
    'nailfold_mask',            true(size(im)),...
    'save_path',                'C:\isbe\nailfold\data\rsa_study\cxx2\vessel_detections\enlargedapex0140_vessel_vd.mat',...
    'plot',                     1);
%%
xi = rand(10,1);
yi = rand(10,1);
num_pts = 10;
sum_x = 0;
sum_x2 = 0;
sum_y = 0;
sum_y2 = 0;

for ii = 1:num_pts

    sum_x = sum_x + xi(ii);
    sum_x2 = sum_x2 + xi(ii)^2;
    sum_y = sum_y + yi(ii);
    sum_y2 = sum_y2 + yi(ii)^2;
end
sigma_x = 1.06*(num_pts^ -.2)*sqrt( (sum_x2 - (sum_x*sum_x/num_pts)) / (num_pts - 1) ); 
sigma_y = 1.06*(num_pts^ -.2)*sqrt( (sum_y2 - (sum_y*sum_y/num_pts)) / (num_pts - 1) );

sigma_x2 = 1.1236*(num_pts^ -.4)*( (sum_x2 - (sum_x*sum_x/num_pts)) / (num_pts - 1) ); 

xm = mean(xi);
sigma_x2 = 1.06*(num_pts^-.2)*sqrt(sum((xi - xm).^2) / (num_pts - 1));
ym = mean(yi);
sigma_y2 = 1.06*(num_pts^-.2)*sqrt(sum((yi - ym).^2) / (num_pts - 1));
%%
xx = repmat(0:31, 32, 1);
yy = xx';
[location_distribution] = build_2d_kernel_distribution(...
           (vh(:,1:2))/8,...
           [xx(:) yy(:)],...
           1-vh(:,4), 0);
 
location_distribution.D_f = reshape(location_distribution.D_f, 32, 32);

%Compute the max y-coordinate for each x-coordinate of the
%candidates, and the displacement to this value
[~, y_max] = max(location_distribution.D_f);
x_i = 8*xx(1,:);
y_i = 8*yy(y_max,1);

candidate_polyfit = interp1(x_i, y_i, vh(:,1), 'linear');
candidate_displacements = vh(:,2) - candidate_polyfit;
%%
%%
load('C:\isbe\nailfold\data\rsa_study\data_lists\image_id_data.mat');
load('C:\isbe\nailfold\data\rsa_study\data_lists\miccai_lists.mat', 'miccai_selection');

im_names = image_id_data.im_names(miccai_selection.validation);
rf = u_load('C:\isbe\nailfold\models\apex\rescoring\corrected_ims_1_228\rf.mat');
predict_apex_candidate_rescore(...
    'image_names',          im_names(229:456),...
    'apex_class_rf',        rf,...
    'data_dir',             [nailfoldroot 'data/rsa_study/master_set/'],...
    'candidates_dir',       'apex_maps\set12g_half_296655\miccai_maxima',...
    'rescore_dir',          'apex_maps\set12g_half_296655\miccai_maxima\rescores_corrected',...
    'hog_dir',              'apex_maps\set12g_half_296655\miccai_maxima\hogs_corrected');

rf = u_load('C:\isbe\nailfold\models\apex\rescoring\corrected_ims_229_456\rf.mat');
predict_apex_candidate_rescore(...
    'image_names',          im_names(1:228),...
    'apex_class_rf',        rf,...
    'data_dir',             [nailfoldroot 'data/rsa_study/master_set/'],...
    'candidates_dir',       'apex_maps\set12g_half_296655\miccai_maxima',...
    'rescore_dir',          'apex_maps\set12g_half_296655\miccai_maxima\rescores_corrected',...
    'hog_dir',              'apex_maps\set12g_half_296655\miccai_maxima\hogs_corrected');
%
compute_apex_candidate_displacements(...
    'image_names', im_names,...
    'data_dir', [nailfoldroot 'data/rsa_study/master_set/'],...
    'vessel_centre_dir', 'vessel_centres\full_centres',...
    'displacement_dir', 'apex_maps\set12g_half_296655\miccai_maxima\displacements_corrected',...
    'candidates_dir', 'apex_maps\set12g_half_296655\miccai_maxima\rescores_corrected',...
    'min_candidates', 3,...
    'initial_thresh', 0.3,...
    'grid_spacing', 8);
%%
load('C:\isbe\nailfold\data\rsa_study\data_lists\image_id_data.mat');
load('C:\isbe\nailfold\data\rsa_study\data_lists\miccai_lists.mat', 'miccai_selection');

make_apex_candidate_class_probs(...
    'image_names',          image_id_data.im_names(miccai_selection.validation),...
    'model_id',             'miccai_class_MAP_corrected',...
    'model_root',           [nailfoldroot,'models/apex'], ...
    'model_name',           'class_map',...
    'data_dir',             [nailfoldroot 'data/rsa_study/master_set/'],...
    'displacement_dir',     'apex_maps\set12g_half_296655\miccai_maxima\displacements_corrected',...
    'candidates_dir',       'apex_maps\set12g_half_296655\miccai_maxima\rescores_corrected',...
    'label_dir',            'apex_maps\set12g_half_296655\miccai_maxima\labels',...
    'plot',                 1,...
    'fig_dir',              []);
%%
load('C:\isbe\nailfold\models\apex\final_MAP\miccai_class_MAP_corrected\class_map.mat');                 
write_array_txt(class_map.x, 'C:\isbe\nailfold\models\apex\final_MAP\miccai_class_MAP_corrected\cxx\class_map_x.txt');
write_array_txt(class_map.y, 'C:\isbe\nailfold\models\apex\final_MAP\miccai_class_MAP_corrected\cxx\class_map_y.txt');
for i_class = 1:3
    write_array_txt(class_map.joint_conditional_probs(:,:,i_class),...
        ['C:\isbe\nailfold\models\apex\final_MAP\miccai_class_MAP_corrected\cxx\class_map_jcb' num2str(i_class) '.txt']);
    write_array_txt(class_map.conditional_class_probs(:,:,i_class),...
        ['C:\isbe\nailfold\models\apex\final_MAP\miccai_class_MAP_corrected\cxx\class_map_ccb' num2str(i_class) '.txt']);
    write_array_txt(class_map.joint_conditional_probs(:,:,i_class),...
        ['C:\isbe\nailfold\models\apex\final_MAP\miccai_class_MAP_corrected\cxx\class_map_pcb' num2str(i_class) '.txt']);
end
write_array_txt(class_map.post_class, 'C:\isbe\nailfold\models\apex\final_MAP\miccai_class_MAP_corrected\cxx\class_map_pc.txt');
write_array_txt(class_map.post_probs, 'C:\isbe\nailfold\models\apex\final_MAP\miccai_class_MAP_corrected\cxx\class_map_pp.txt');