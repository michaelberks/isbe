%--------------------------------------------------------------------------
% ************* Script for Chitra - UG project option 2011 ****************
%--------------------------------------------------------------------------
%% 1. Load in density data and get list of mammograms to use

%Get data from Camilla's excel spreadsheet (2009 UG projects)
[lcc_nums] = xlsread('C:\isbe\ug_project_options\2011\chitra\area gland 5% compilation.xls', 1, 'A2:A169');
[lcc_dens] = xlsread('C:\isbe\ug_project_options\2011\chitra\area gland 5% compilation.xls', 1, 'J2:J169');

%Discrad 16th entry as no reading from Dr Wilson
lcc_nums(16) = [];
lcc_dens(16) = [];

%Discard entries with mammo names over 1000
discard_idx = lcc_nums > 999;
lcc_nums(discard_idx) = [];
lcc_dens(discard_idx) = [];

%Get list of mammo names
num_mammos = length(lcc_nums);
mam_names = cell(num_mammos,1);
for ii = 1:num_mammos
    mam_names{ii} = [zerostr(lcc_nums(ii),3) 'LCC'];
end

%Get list of matching mammograms and segmentations
[tif_names missing_tif] = match_mammo_names('I:\Density_feasibility_study_HD\', mam_names, [], [], '.tif');
[seg_names missing_seg] = match_mammo_names('J:\stepwedge\segmentation\', mam_names, [], [], '.mat');
%--------------------------------------------------------------------------
%% 2. Process mammos separately first

%Make new directory
mkdir C:\isbe\ug_project_options\2011\chitra\small_mammos;
mkdir C:\isbe\ug_project_options\2011\chitra\segmentations;

%Loop through mammos on hard drive reducing bit depth and resizing
for ii = 1:num_mammos
    mammo = imread(['I:\density_feasibility_study_HD\' tif_names{ii}]);
    mammo = uint8(mammo / 256);
    
    %Check mammogram is portrait orientated
    [r c] = size(mammo);
    if c > r
        mammo = rot90(mammo);
    end

    % resize to have 1024 rows    
    mammo = imresize(mammo, [1024 nan], 'nearest');
    [r c] = size(mammo);
    
    %check if right or left
    right = mam_names{ii}(4) == 'R';
    
    %check orientation of mammogram
    r1 = round(r / 3);
    r2 = 2*r1;
    c_half = round(c / 2);
    
    mam_mean = mean(mammo(:));
    centre_mean = mean(mammo(r1:r2,:));
    
    sum1 = sum(centre_mean(1:c_half) > mam_mean);
    sum2 = sum(centre_mean(c_half+1:end) > mam_mean);
    
    %if sum1 > sum2 then the breast is on the left of the image
    if (sum1 > sum2) == right
        mammo = rot90(mammo, 2);
    end
    
    save(['C:\isbe\ug_project_options\2011\chitra\small_mammos\' mam_names{ii} '.mat'], 'mammo');
    clear mammo;

    %Copy over the segmentation too
    copyfile(...
        ['J:\stepwedge\segmentation\' seg_names{ii}],...
        ['C:\isbe\ug_project_options\2011\chitra\segmentations\' seg_names{ii}]);
end
%--------------------------------------------------------------------------
%% 3. Create multiblob synthetic mammograms
%--------------------------------------------------------------------------

%pre-allocate storage for densities
density_a = zeros(num_mammos, 1);
density_b = zeros(num_mammos, 1);

%Define grey levels for fat and gland
fat_g = 120;
gland_g = 196;
    
%define paths to load/save data
blob_path = 'C:\isbe\ug_project_options\2011\chitra\blobs\blob';
grey_path = 'C:\isbe\ug_project_options\2011\chitra\blobs_grey\blob';
mammo_path = 'C:\isbe\ug_project_options\2011\chitra\small_mammos\';
seg_path = 'C:\isbe\ug_project_options\2011\chitra\segmentations\';

%Load mass model
load D:\isbe\dev\mass_model\models\model_w500_50K.mat

%prepare fields from mass model
%--------------------------------------------------------------------------
mean_shape  = mass_model.mean_shape;
P_shape     = mass_model.P_shape;
L_shape     = mass_model.L_shape;
mean_scale  = mass_model.mean_scale;
P_scale     = mass_model.P_scale;
mean_tex    = mass_model.mean_tex;
P_tex       = mass_model.P_tex;
L_tex       = mass_model.L_tex;

P_com         = mass_model.P_com;
L_com         = mass_model.L_com;

W_shape     = mass_model.W_shape;
W_tex       = mass_model.W_tex;
W_scale     = mass_model.W_scale;

mean_shape_pl = mass_model.mean_shape_pl;
size_shape_vec = length(mean_shape) / 2;

k_shape     = length(L_shape);
k_tex       = length(L_tex);

scaled_mean_shape = mass_model.scaled_mean_shape;
  
%Define source points for TPS - as row vectors
s_x = scaled_mean_shape(1:size_shape_vec);
s_y = scaled_mean_shape(size_shape_vec+1:end);

%Define points to be interpolated by TPS - as row vectors
i_x = mean_shape_pl(:,1)';
i_y = mean_shape_pl(:,2)';

len = 500;
keep = 10;
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
mam_count = 1;
for dd = 1:num_mammos
    display(['**************** Blob ' zerostr(dd,3) ' *********************']);
    
    %Randomly sample a density change
    density_change = (randn * 10 / 2) + 10;
    
    %Check the mammogram is dense enough to begin with:
    if density_change > lcc_dens(dd)
        display(['Skipping bob ' zerostr(dd,3)]);
        continue;
    end
    %Load in mammogram and segmentation
    mammogram = u_load([mammo_path mam_names{dd}]);
    segmentation = u_load([seg_path seg_names{dd}]);
    mammogram = imresize(mammogram, segmentation.size, 'nearest');

    %Get breast shape as continous polygon
    dists = cumsum([0; sqrt(sum(diff(segmentation.breast_border).^2, 2))]);
    breast_shape = interp1(dists, segmentation.breast_border, linspace(0, dists(end), len), 'linear');
    
    %Make mask of breast and compute breast area
    breast_mask = roipoly(mammogram, breast_shape(:,1), breast_shape(:,2));
    breast_area = polyarea(breast_shape(:,1), breast_shape(:,2));
    
    %Workout furthest point to right and assume nipple position to locate
    %first blob
    [dummy nipple_idx] = max(breast_shape(:,1));
    start_point = breast_shape(nipple_idx,:) - [25 0];

    %Workout initial blob shape to match breast shape
    [d, aligned_shape, transform] = mb_procrustes(mass_model.mean_target,breast_shape);
    b_shape = P_shape' * (aligned_shape(:)' - mean_shape)';
    

    %Check if density change is positive or negative, this will detrmine
    %whether the first mammogram is A or B
    if density_change > 0
        blob_name_1 = [blob_path zerostr(mam_count,3) '_a'];
        blob_name_2 = [blob_path zerostr(mam_count,3) '_b'];
        
        a_idx = 1;
        b_idx = 2;
        
        blob_sum = breast_area * lcc_dens(dd) / 100;
    else
        blob_name_1 = [blob_path zerostr(mam_count,3) '_b'];
        blob_name_2 = [blob_path zerostr(mam_count,3) '_a'];
        
        a_idx = 2;
        b_idx = 1;
        
        blob_sum = breast_area * (lcc_dens(dd)-density_change) / 100;
    end
    
    %Pre-allocate blob mask and blob mammogram
    blob_mask = false(size(breast_mask));
    blob_mammogram = zeros(size(breast_mask));
    num_blobs = 0;
    
    %Now sample blobs to create first blob mammogram
    while sum(blob_mask(:)) < blob_sum
        
        %Increment blob count
        num_blobs = num_blobs + 1;
        
        %Generate random rotation of shape
        theta = 2*pi*rand;
        theta_m = [cos(theta) sin(theta); -sin(theta) cos(theta)];
        
        % Compute new shape vector, texture vector and scale
        blob_com = (randn(length(L_com), 1) .* sqrt(L_com));

        Q_shape = P_com(1:k_shape,:); 
        Q_tex = P_com(k_shape+1:k_shape + k_tex,:);
        Q_scale = P_com(end, :);

        %Use sampled parameters to generate texture
        blob_tex = mean_tex + (P_tex*inv(W_tex)*Q_tex*blob_com)';
        blob_tex = blob_tex * 1.2;
        
        %For first blob, we use the initial shape to match the breast shape
        %and locate just inward of the nipple
        if num_blobs == 1
            b_shape(keep+1:end) = (randn(length(mass_model.L_shape)-keep, 1) .* sqrt(mass_model.L_shape(keep+1:end)));

            blob_shape = reshape(mass_model.P_shape*b_shape + mass_model.mean_shape', [], 2);

            blob_area = polyarea(blob_shape(:,1), blob_shape(:,2));

            scale_factor = sqrt(blob_sum /(blob_area*2));

            blob_shape = scale_factor * blob_shape* inv(transform.T);
            blob_shape(:,1) = blob_shape(:,1) - blob_shape(nipple_idx,1) + start_point(1);
            blob_shape(:,2) = blob_shape(:,2) - blob_shape(nipple_idx,2) + start_point(2);
            
        else
        
            %Use sampled parameters to generate shape
            blob_shape = mean_shape + (P_shape*inv(W_shape)*Q_shape*blob_com)';
            blob_shape = rand*reshape(blob_shape, [], 2)* theta_m;

            %Workout where we can stick this blob
            blob_size = round(max(max(blob_shape) - min(blob_shape)));
            blob_map = imerode(breast_mask, strel('disk', blob_size)) & ~blob_mask;
            [pts_y pts_x] = find(blob_map);

            %Are there any viable pts...?
            if isempty(pts_x); continue; end
            
            %If yes, select one at random
            centre_idx = ceil(length(pts_x)*rand);

            %Centre the new blow shape here
            blob_shape(:,1) = blob_shape(:,1) - mean(blob_shape(:,1)) + pts_x(centre_idx);
            blob_shape(:,2) = blob_shape(:,2) - mean(blob_shape(:,2)) + pts_y(centre_idx);
        end
        
        %Get a mask of the new blob, and compute the indices of pixels
        %belonging to the blob
        blob_bw = roipoly(blob_mask, blob_shape(:,1), blob_shape(:,2));
        blob_idx = find(blob_bw);
        [blob_y blob_x] = ind2sub(size(blob_mask), blob_idx);

        % Compute TPS warp to texture from mean to new blob
        %%%%

        %Define displacement to target points
        z_x = blob_shape(:,1)';
        z_y = blob_shape(:,2)';

        %Compute displacement of interpolated points        
        T = geom_alignpoints_2d([z_x; z_y], [s_x; s_y], [],...
            'transform', 'spline');
        [f_xy] = geom_transformpoints([i_x; i_y], T);

        %Interpolate uneven grid of sampled texture points to pixel
        %lattice
        blob_shape_tex = griddata(f_xy(1,:), f_xy(2,:), blob_tex,...
            blob_x, blob_y);

        %Remove any unrealistic values
        blob_shape_tex(isnan(blob_shape_tex)) = 0;
        blob_shape_tex(blob_shape_tex < 0) = 0;

        %Add the new blob to the blob mask and blob mammogram
        blob_mammogram(blob_bw) =  blob_mammogram(blob_bw)+blob_shape_tex;
        blob_mask = blob_mammogram > 0;
       
    end    
    
    %Save the grey-scale synthetic mammogram
    save([grey_path zerostr(mam_count,3) '.mat'], 'blob_mammogram');
    
    %Workout density of initial blob mammo
    densities(1) = sum(blob_mask(:)) / breast_area;
    blob_sum = (densities(1) - abs(density_change)/100)*breast_area;
    
    %Create blob base base
    blob_base = mammogram;
    blob_base(breast_mask) = fat_g;
    
    %create first blob mammogram
    blob_mammo_1 = blob_base;
    blob_mammo_1(blob_mask) = gland_g;
    
    %Now loop through
    thresh = 0;
    thresh_mask = blob_mammogram > 0;
    while sum(thresh_mask(:)) >= blob_sum
        thresh = thresh + 1;
        thresh_mask = blob_mammogram > thresh;
    end
    
    %Workout density of second blob mammo
    densities(2) = sum(thresh_mask(:)) / breast_area;
    
    %create second blob mammogram
    blob_mammo_2 = blob_base;
    blob_mammo_2(thresh_mask) = gland_g;
    
    %Save the 2 densities
    density_a(dd) = densities(a_idx);
    density_b(dd) = densities(b_idx);

    %Save the 2 blob mammograms as jpegs/mat files
    write_im_from_colormap(blob_mammo_1, [blob_name_1 '.jpg'], gray(256));
    write_im_from_colormap(blob_mammo_2, [blob_name_2 '.jpg'], gray(256));
    save([blob_name_1 '.mat'], 'blob_mammo_1');
    save([blob_name_2 '.mat'], 'blob_mammo_2');
    
    if dd < 11
        figure;
        subplot(1,2,1); imagesc(blob_mammo_1); axis image; colormap(gray(256));
        title(['Density = ' num2str(densities(1))]);
        subplot(1,2,2); imagesc(blob_mammo_2); axis image; colormap(gray(256));
        title(['Density = ' num2str(densities(2))]);
    end
    mam_count = mam_count + 1;
    
end
save C:\isbe\ug_project_options\2011\chitra\blob_densities density_a density_b
%%
densities = zeros(50,2);
fat_g = 120;
gland_g = 196;
ab = 'ab';
for ii = 1:50
    for kk = 1:2
        blob = u_load(['C:\isbe\ug_project_options\2011\chitra\blobs\blob' zerostr(ii,3) '_' ab(kk) '.mat']);
        fat_mask = blob == fat_g;
        gland_mask = blob == gland_g;
        mam_label = bwlabel(fat_mask | gland_mask,4);
        num_labels = max(mam_label(:));
        counts = zeros(num_labels,1);
        for jj = 1:num_labels
            counts(jj) = sum(mam_label(:) == jj);
        end
        [d max_l] = max(counts);
        mam_mask = mam_label == max_l;
        
        gland_mask = gland_mask & mam_mask;
        
        %figure; imagesc(gland_mask); axis image;
        densities(ii,kk) = 100 * sum(gland_mask(:)) / sum(mam_mask(:));
    end
end
xlswrite('temp.xls', densities);    
%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%4. Make VAS forms
for reader_id = 1:15
    form_folder = ['C:\isbe\ug_project_options\2011\chitra\forms\reader' zerostr(reader_id,2) '\'];
    mkdir(form_folder);
    
    for form_type = 1:4
        if form_type == 1
            end_idx = 20;
        else
            end_idx = 10;
        end
        for form_idx = 1:end_idx
            make_form_pdf_chitra(form_type, form_idx, reader_id, form_folder);
        end
    end
end
%%
for reader_id = 1:15
    form_folder = ['C:\isbe\ug_project_options\2011\chitra\forms\reader' zerostr(reader_id,2) '\'];
    mkdir(form_folder);
    
    for form_idx = 1:20
        make_form_pdf_chitra(1, form_idx, reader_id, form_folder);
    end
    for form_idx = 1:10
        make_form_pdf_chitra2(form_idx, reader_id, form_folder);
    end
end
%%
pdf_list = dir('C:\isbe\ug_project_options\2011\chitra\scanned\*.pdf');
pdf_dir = 'C:\isbe\ug_project_options\2011\chitra\scanned\';
xls_file = 'C:\isbe\ug_project_options\2011\chitra\results\results.xls';
for ii = 1:length(pdf_list)
    read_chitra_form_batch([pdf_dir pdf_list(ii).name], [], xls_file, 0);
end
%%
xls_file = 'C:\isbe\ug_project_options\2011\chitra\results\results.xls';
user_path = 'C:\isbe\ug_project_options\2011\chitra\users\user';
for ii = 10
    if exist([user_path zerostr(ii,2) '_1.mat'], 'file')
        user_data = u_load([user_path zerostr(ii, 2) '_1.mat']);
        xlswrite(xls_file, user_data.blob_order, ii, 'A2:A101');
        xlswrite(xls_file, user_data.blob_timings, ii, 'C2:C101');
    end
    if exist([user_path zerostr(ii,2) '_2.mat'], 'file')
        user_data = u_load([user_path zerostr(ii, 2) '_2.mat']);
        xlswrite(xls_file, user_data.blob_order, ii, 'D2:D51');
        xlswrite(xls_file, user_data.blob_timings, ii, 'G2:G51');
    end
%     xlswrite(xls_file, {'Time'}, ii, 'C1');
%     xlswrite(xls_file, {'Time'}, ii, 'G1');
%     xlswrite(xls_file, {'Time'}, ii, 'J1');
%     xlswrite(xls_file, {'Image order'}, ii, 'H1');
%     xlswrite(xls_file, {'Image order'}, ii, 'H1');
%     xlswrite(xls_file, {'Density change'}, ii, 'I1');
    if exist([user_path zerostr(ii,2) '_3.mat'], 'file')
        user_data = u_load([user_path zerostr(ii, 2) '_3.mat']);
        xlswrite(xls_file, user_data.blob_order, ii, 'H2:H51');
        xlswrite(xls_file, user_data.ratings, ii, 'I2:I51');
        xlswrite(xls_file, user_data.blob_timings, ii, 'J2:J51');
    end
end
%%
fid = fopen('C:\isbe\ug_project_options\2011\chitra\software\missing_images.txt');
def_txt = textscan(fid, '%s %s', 'delimiter', '=');
fclose(fid);

for ii = 1:2
	missing_ims = ...
        str2num(def_txt{2}{strncmpi(def_txt{1}, ['reader_' zerostr(ii,2)], 9)});
	display(missing_ims);
end
    def_txt{2}{strncmpi(def_txt{1}, 'xls_names_dir', length('xls_names_dir'))};