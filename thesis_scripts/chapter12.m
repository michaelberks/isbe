% Chapter 12, a few results for future work

%Build model data for set of unused masses...
load C:\isbe\dev\files\u_files.mat
load C:\isbe\dev\mass_model\models\model_w500_50K.mat

unused_idx = setdiff(1:179, idx_u1);
mass_list = dir('C:\isbe\dev\masses\*.mat');
unused_list = mass_list(unused_idx);
[unused_model] = get_unused_mass_parameters('C:\isbe\dev\mass_model\models\model_w500_50K.mat', unused_list);
save C:\isbe\dev\mass_model\models\unused_masses_model unused_model
%%
pairs_idx = [];
for ii = 1:length(unused_list)
    mass_stem = unused_list(ii).name(5:end-4);
    if ~isempty(strfind(mass_stem, 'CC'))
        mass_stem(5:6) = 'ML';
    else
        mass_stem(5:6) = 'CC';
    end
    for jj = 1:length(u_files1);
        mass_stem2 = u_files1(jj).name(5:end-4);
        if strcmp(mass_stem, mass_stem2)
            display([unused_list(ii).name, ' ', u_files1(jj).name]);
            pairs_idx = [pairs_idx; ii jj]; %#ok
        end
    end
end
%%
figure; hold on;
for ii = 1:size(pairs_idx,1);
    x = [unused_model.B_com(1,pairs_idx(ii,1)) mass_model.B_com(1,pairs_idx(ii,2))];
    y = [unused_model.B_com(2,pairs_idx(ii,1)) mass_model.B_com(2,pairs_idx(ii,2))];
    %plot(x,y);
    
    z = [unused_model.B_com(3,pairs_idx(ii,1)) mass_model.B_com(3,pairs_idx(ii,2))];
    plot3(x,y,z);
end
%%
N = size(mass_model.B_com,2);
Nu = size(pairs_idx,1);
dist = zeros(Nu,1);
avg_dist = zeros(Nu,1);

for ii = 1:Nu;
    b1 = unused_model.B_com(:,pairs_idx(ii,1));
    b2 = mass_model.B_com(:,pairs_idx(ii,2));
    
    dist(ii) = sqrt(sum((b1-b2).^2));
    avg_dist(ii) = mean(sqrt(sum((repmat(b1,1,N) - mass_model.B_com).^2,2)));
    display(['View dist = ', num2str(dist(ii)), ' Avg dist = ', num2str(avg_dist(ii))]);
end
%%
for ii = 1:10%size(pairs_idx,1);
    mass1 = u_load(['C:\isbe\dev\masses\', unused_list(pairs_idx(ii,1)).name]);
    mass2 = u_load(['C:\isbe\dev\masses\', u_files1(pairs_idx(ii,2)).name]);
    
    figure;
    subplot(1,2,1); imagesc(mass1.subtract_ROI); colormap(gray(256)); axis image; hold on;
    plot(mass1.mass_outline(:,1), mass1.mass_outline(:,2));
    subplot(1,2,2); imagesc(mass2.subtract_ROI); colormap(gray(256)); axis image; hold on;
    plot(mass2.mass_outline(:,1), mass2.mass_outline(:,2));
    
    clear mass1 mass2
end
%%
b_views = [mass_model.B_com(:,pairs_idx(:,2)); unused_model.B_com(:,pairs_idx(:,1))];
k_view = size(b_views,1);
cov_views = cov(b_views');
for ii = 1:10
    b1 = unused_model.B_com(:,pairs_idx(ii,1));
    b2 = mass_model.B_com(:,pairs_idx(ii,2));
    b_view = nan(1,k_view);
    b_view(1:end/2) = b2;
    [c_mean c_covar] = condition_gaussian(zeros(1,k_view), cov_views, b_view);
    
    figure; hold on;
    plot(b2(1), b2(2), 'b+');
    plot(b1(1), b1(2), 'gx');
    x_axis = c_covar(1,1:2) / sqrt(sum(c_covar(1,1:2).^2));
    r_x = 2*sqrt(c_covar(1,1));
    r_y = 2*sqrt(c_covar(2,2));
    plot_ellipse(r_x, r_y, x_axis, c_mean(1), c_mean(2));
end
%%
subtract_mass_it(mass_files(1:10), 'C:\isbe\dev\annotations\', 100, 20, 10, 18, 1, 'biharmTPS', 1);
%%
for ii = 1:length(unused_list)
    mass_stem = unused_list(ii).name(5:end-4);
    if ~isempty(strfind(mass_stem, 'CC'))
        mass_stem(5:6) = 'ML';
        cc_list(ii).name = unused_list(ii).name;
        
        for jj = 1:length(u_files1);
            mass_stem2 = u_files1(jj).name(5:end-4);
            if strcmp(mass_stem, mass_stem2)
                mlo_list(ii).name = u_files1(jj).name;
                display([cc_list(ii).name, ' ', mlo_list(ii).name]);
            end
        end
    else
        mass_stem(5:6) = 'CC';
        mlo_list(ii).name = unused_list(ii).name;
        
        for jj = 1:length(u_files1);
            mass_stem2 = u_files1(jj).name(5:end-4);
            if strcmp(mass_stem, mass_stem2)
                cc_list(ii).name = u_files1(jj).name;
                display([cc_list(ii).name, ' ', mlo_list(ii).name]);
            end
        end
    end
    
end
%%
[cc_model cc_id]...
    = generate_mass_AM(cc_list, 'C:\isbe\dev\mass_model\models\paired_cc_model');
[mlo_model mlo_id]...
    = generate_mass_AM(mlo_list, 'C:\isbe\dev\mass_model\models\paired_mlo_model');