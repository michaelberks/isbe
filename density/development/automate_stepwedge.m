%Script looking at methods for automating the location of the stepwedge
%from the training set of mammograms

%First step, let's create a dataset of image regions containing the wedge
mam_list = dir('G:\*1824*.tif');
mkdir C:\isbe\density\stepwedges
resize_factor = 0.176; %as used in main program
%
for ii = 1:30;
    
    small_size = isempty(strfind(mam_list(ii).name, '2430'));
    mam = imresize(imread(['G:\' mam_list(ii).name]), resize_factor);
    mam = medfilt2(mam);
    mam = 2^16 - mam - 1;
    mam = double(mam) ./ max(double(mam(:)));
    if small_size
        mam = rot90(mam);
        y_start = 150;
        y_end = 200;
    else %filmsize == 2 %24x30 mammogram
        y_start = 250;
        y_end = 400;      
    end
    
    pos = [1 y_start 400 size(mam,1)-y_end];
    step_image = mam(pos(2):pos(2)+pos(4)-1, pos(1):pos(1)+pos(3)-1);
    imwrite(step_image, ['C:\isbe\density\stepwedges\' mam_list(ii).name(1:end-3) 'bmp']);
end
%%
step_list = dir('C:\isbe\density\stepwedges\images\*.bmp');
for ii = 1:length(step_list);
    display(step_list(ii).name);
end
%%
%Compute points for each stepwedge image using model (either 1 or 2,
%depending on params file). The params file contains the list of images
%(containing stepwedges) to process
cmd = 'C:\isbe\cxx\vxl\bin\isbe_autoland\urf\tools\release\urf_graph_model_search.exe';
options = ' -p C:\isbe\density\stepwedges\params\urf_graph_model_search_params2.txt';
[status] = system([cmd options]);
%%
%Have a look at the points predicted by the model on some images
for ii = [14 15 56 58]%1:20%length(step_list);
    
    %load in the stepwedge region
    step_im = imread(['C:\isbe\density\stepwedges\images\', step_list(ii).name]);
    figure;
    
    %load in the stepwedge points computed by model 1
    fid = fopen(['C:\isbe\density\stepwedges\points\model1\', step_list(ii).name, '.pts']);
    step_pts_txt = textscan(fid, '%s %s');
    fclose(fid);
    x_pts = str2double(step_pts_txt{1}(4:16));
    y_pts = str2double(step_pts_txt{2}(4:16));
    
    %Plot these points on an image of the stepwedge
    subplot(1,2,1); imagesc(step_im); axis image; colormap(gray(256)); hold on;
    plot(x_pts, y_pts, 'gx');
    
    %load in the stepwedge points computed by model 2
    fid = fopen(['C:\isbe\density\stepwedges\points\model2\', step_list(ii).name, '.pts']);
    step_pts_txt = textscan(fid, '%s %s');
    fclose(fid);
    x_pts = str2double(step_pts_txt{1}(4:16));
    y_pts = str2double(step_pts_txt{2}(4:16)); 
    
    %plot these points on the image
    subplot(1,2,2); imagesc(step_im); axis image; colormap(gray(256)); hold on;
    plot(x_pts, y_pts, 'rx');
%     plot(x_pts(5:6), y_pts(5:6), 'r');
%     plot(x_pts(7:8), y_pts(7:8), 'b');
%     plot(x_pts(9:10), y_pts(9:10), 'g');
%     plot(x_pts(11:12), y_pts(11:12), 'y');
end

%%
ii = 58;
    
%load in the stepwedge region
step_im = imread(['C:\isbe\density\stepwedges\images\', step_list(ii).name]);
f1 = figure(...
    'windowstyle', 'normal',...
    'Units', 'centimeters ',...
    'position', [0 0 10.0 20.0],...
    'PaperPositionMode','auto');

axes(...
    'Units', 'centimeters ',...
    'position', [0 0 10.0 20.0]);

imagesc(step_im); axis image off; colormap(gray(256)); hold on;

%load in the stepwedge points computed by model 1
fid = fopen(['C:\isbe\density\stepwedges\points\model1\', step_list(ii).name, '.pts']);
step_pts_txt = textscan(fid, '%s %s');
fclose(fid);
x_pts = str2double(step_pts_txt{1}(4:16));
y_pts = str2double(step_pts_txt{2}(4:16));

%Plot these points on an image of the stepwedge 
plot(x_pts, y_pts, 'gx', 'markersize', 8);

%load in the stepwedge points computed by model 2
fid = fopen(['C:\isbe\density\stepwedges\points\model2\', step_list(ii).name, '.pts']);
step_pts_txt = textscan(fid, '%s %s');
fclose(fid);
x_pts = str2double(step_pts_txt{1}(4:16));
y_pts = str2double(step_pts_txt{2}(4:16)); 

%plot these points on the image
plot(x_pts, y_pts, 'ro', 'markersize', 8);
print('-dtiff', '-noui', '-painters', f1, '-r300', ['C:\isbe\density\stepwedges\figures\sw' zerostr(ii,2) '_models1+2.tif']);

%%
%Align the stepwedge points predicted by the models for each stepwedge, and
%compute a mean shape using Procrustes alignment
valid_steps = 1:58;%setdiff(1:length(step_list), [15 56 58]);

step_shapes = zeros(length(valid_steps),52);
step_areas = zeros(length(valid_steps),1);

%Load in points for each stepwedge and store in shape matrix
for jj = 1:length(valid_steps);
    ii = valid_steps(jj);
    fid = fopen(['C:\isbe\density\stepwedges\points\model1\', step_list(ii).name, '.pts']);
    step_pts_txt = textscan(fid, '%s %s');
    fclose(fid);
    
    x_pts1 = str2double(step_pts_txt{1}(4:16));
    y_pts1 = str2double(step_pts_txt{2}(4:16));
    
    fid = fopen(['C:\isbe\density\stepwedges\points\model2\', step_list(ii).name, '.pts']);
    step_pts_txt = textscan(fid, '%s %s');
    fclose(fid);
    x_pts2 = str2double(step_pts_txt{1}(4:16));
    y_pts2 = str2double(step_pts_txt{2}(4:16));
    
     
    step_shapes(jj, :) = [x_pts1' x_pts2' y_pts1' y_pts2'];
    step_areas(jj) = polyarea([x_pts1' x_pts2'], [y_pts1' y_pts2']);    

end
%Use procrustes analysis to align shapes
args.seed = 1;
args.shiftOrigin = 0;
args.area = mean(step_areas);
[a_shapes, a_scales, step_mean, a_rots, a_trans a_origins]...
    = align_shapes(step_shapes, args);

% For each stepwedge compute the alignment error between itself and the
% mean shape
align_error = zeros(58,1);
for ii = 1:58;
    fid = fopen(['C:\isbe\density\stepwedges\points\model1\', step_list(ii).name, '.pts']);
    step_pts_txt = textscan(fid, '%s %s');
    fclose(fid);
    
    x_pts1 = str2double(step_pts_txt{1}(4:16));
    y_pts1 = str2double(step_pts_txt{2}(4:16));
    
    fid = fopen(['C:\isbe\density\stepwedges\points\model2\', step_list(ii).name, '.pts']);
    step_pts_txt = textscan(fid, '%s %s');
    fclose(fid);
    x_pts2 = str2double(step_pts_txt{1}(4:16));
    y_pts2 = str2double(step_pts_txt{2}(4:16));
    
    [align_error(ii),Z,transform] = mb_procrustes(step_mean, [x_pts1 y_pts1; x_pts2 y_pts2]);
end

%%
for ii = [14 15 56 58];
    fid = fopen(['C:\isbe\density\stepwedges\points\model1\', step_list(ii).name, '.pts']);
    step_pts_txt = textscan(fid, '%s %s');
    fclose(fid);
    
    x_pts1 = str2double(step_pts_txt{1}(4:16));
    y_pts1 = str2double(step_pts_txt{2}(4:16));
    
    fid = fopen(['C:\isbe\density\stepwedges\points\model2\', step_list(ii).name, '.pts']);
    step_pts_txt = textscan(fid, '%s %s');
    fclose(fid);
    x_pts2 = str2double(step_pts_txt{1}(4:16));
    y_pts2 = str2double(step_pts_txt{2}(4:16));

    step_shape = [x_pts1 y_pts1; x_pts2 y_pts2]; 
    [error1,step_mean_a] = mb_procrustes(step_shape, step_mean);
    dists = sum((step_mean_a - step_shape).^2,2);
    [dummy gd_pts] = sort(dists, 'ascend');
    
    n = 13;
    [error2,step_mean_a2, transform] = mb_procrustes(step_shape(gd_pts(1:n),:), step_mean(gd_pts(1:n),:));
    
    step_mean_a3 = transform.b * step_mean * transform.T + repmat(transform.c(1,:),26,1);
    
    step_im = imread(['C:\isbe\density\stepwedges\images\', step_list(ii).name]);
    figure;    
    subplot(1,2,1); imagesc(step_im); axis image; colormap(gray(256)); hold on;
    plot(step_shape(1:13,1), step_shape(1:13,2), 'rx');
    plot(step_shape(14:26,1), step_shape(14:26,2), 'r+');
    plot(step_mean_a(1:13,1), step_mean_a(1:13,2), 'yx');
    plot(step_mean_a(14:26,1), step_mean_a(14:26,2), 'y+');
    for jj = 1:5
        plot(step_shape(gd_pts(1:n),1), step_shape(gd_pts(1:n),2), 'bo');
    end
    
    subplot(1,2,2); imagesc(step_im); axis image; colormap(gray(256)); hold on;
    plot(step_shape(gd_pts(1:n),1), step_shape(gd_pts(1:n),2), 'rx');
    plot(step_mean_a3(:,1), step_mean_a3(:,2), 'y+');
    
    display(['Original error = ', num2str(error1), ' Improved error = ', num2str(error2)]);
end
%%
%For each stepwedge, align the mean stepwedge shape to the predicted
%stepwedge shape and use the aligned mean shape - the theory being this has
%a more robust shape than the one predicted for a single image
for ii = 1:20;
    
    %load in model 1 points
    fid = fopen(['C:\isbe\density\stepwedges\points\model1\', step_list(ii).name, '.pts']);
    step_pts_txt = textscan(fid, '%s %s');
    fclose(fid);
    
    x_pts1 = str2double(step_pts_txt{1}(4:16));
    y_pts1 = str2double(step_pts_txt{2}(4:16));
    
    %load in model 2 points
    fid = fopen(['C:\isbe\density\stepwedges\points\model2\', step_list(ii).name, '.pts']);
    step_pts_txt = textscan(fid, '%s %s');
    fclose(fid);
    x_pts2 = str2double(step_pts_txt{1}(4:16));
    y_pts2 = str2double(step_pts_txt{2}(4:16));

    %concatenate to form predicted model shape the align mean stepwedge
    %shape to this
    step_shape = [x_pts1 y_pts1; x_pts2 y_pts2]; 
    [error1,step_mean_a] = mb_procrustes(step_shape, step_mean);
    
    %display the stepwedge with the model points (red) and aligned mean
    %points (yellow) plotted
    step_im = imread(['C:\isbe\density\stepwedges\images\', step_list(ii).name]);
    figure;    
    imagesc(step_im); axis image; colormap(gray(256)); hold on;
    plot(step_shape(1:13,1), step_shape(1:13,2), 'rx');
    plot(step_shape(14:26,1), step_shape(14:26,2), 'r+');
    plot(step_mean_a(1:13,1), step_mean_a(1:13,2), 'yx');
    plot(step_mean_a(14:26,1), step_mean_a(14:26,2), 'y+');
end
%%
%Use the predicted stepwedge points to predict the centre position of each step
step_num = [24 29 31 35 39]';
step_pts1 =[x_pts1(12)                    y_pts1(12);...
            mean([x_pts1(10) x_pts2(13)]) mean([y_pts1(10) y_pts2(13)]);...
            mean([x_pts1(8) x_pts2(11)])  mean([y_pts1(8) y_pts2(11)]);...
            x_pts2(9)                     y_pts2(9);...
            mean([x_pts1(6) x_pts2(7)])   mean([y_pts1(6) y_pts2(7)])];
step_pts2 =[x_pts1(13)                    y_pts1(13);...
            x_pts1(11)                    y_pts1(11);...
            mean([x_pts1(9) x_pts2(12)])  mean([y_pts1(9) y_pts2(12)]);...
            x_pts2(10)                     y_pts2(10);...
            mean([x_pts1(7) x_pts2(8)])   mean([y_pts1(7) y_pts2(8)])];

step_pts1i = interp1(step_num, step_pts1, 1:39, 'linear', 'extrap');
step_pts2i = interp1(step_num, step_pts2, 1:39, 'linear', 'extrap');

%Display the stepwedge with the predicted steps
figure; imagesc(step_im); axis image; colormap(gray(256)); hold on;
plot([step_pts1i(:,1) step_pts2i(:,1)]', [step_pts1i(:,2) step_pts2i(:,2)]');
%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%as above but for each stepwedge
for ii = [1 2 14]%1:20%[14 15 50:58];
    
    %Load in the points from model 1
    fid = fopen(['C:\isbe\density\stepwedges\points\model1\', step_list(ii).name, '.pts']);
    step_pts_txt = textscan(fid, '%s %s');
    fclose(fid);
    
    x_pts1 = str2double(step_pts_txt{1}(4:16));
    y_pts1 = str2double(step_pts_txt{2}(4:16));
    
    %load in the points from model 2
    fid = fopen(['C:\isbe\density\stepwedges\points\model2\', step_list(ii).name, '.pts']);
    step_pts_txt = textscan(fid, '%s %s');
    fclose(fid);
    x_pts2 = str2double(step_pts_txt{1}(4:16));
    y_pts2 = str2double(step_pts_txt{2}(4:16));

    %Concatenate to form step shape and align mean shape to model shape
    step_shape = [x_pts1 y_pts1; x_pts2 y_pts2]; 
    [align_error,step_mean_a] = mb_procrustes(step_shape, step_mean);
    
    %try and do the good points corrections for shapes with high alignment
    %error
    if align_error > 0.001
        dists = sum((step_mean_a - step_shape).^2,2);
        [dummy gd_pts] = sort(dists, 'ascend');
        
        n = 13;
        [error2,step_mean_a2, transform] = mb_procrustes(step_shape(gd_pts(1:n),:), step_mean(gd_pts(1:n),:));
        step_mean_a3 = transform.b * step_mean * transform.T + repmat(transform.c(1,:),26,1);
        
        x_pts1 = step_mean_a3(1:13,1);
        x_pts2 = step_mean_a3(14:26,1);
        y_pts1 = step_mean_a3(1:13,2);
        y_pts2 = step_mean_a3(14:26,2);
    end
    
    %Predict the poisition of each step centre from the aligned shape
    step_num = [24 29 31 35 39]';
    step_pts1 =[x_pts1(12)                    y_pts1(12);...
                mean([x_pts1(10) x_pts2(13)]) mean([y_pts1(10) y_pts2(13)]);...
                mean([x_pts1(8) x_pts2(11)])  mean([y_pts1(8) y_pts2(11)]);...
                x_pts2(9)                     y_pts2(9);...
                mean([x_pts1(6) x_pts2(7)])   mean([y_pts1(6) y_pts2(7)])];
    step_pts2 =[x_pts1(13)                    y_pts1(13);...
                x_pts1(11)                    y_pts1(11);...
                mean([x_pts1(9) x_pts2(12)])  mean([y_pts1(9) y_pts2(12)]);...
                x_pts2(10)                     y_pts2(10);...
                mean([x_pts1(7) x_pts2(8)])   mean([y_pts1(7) y_pts2(8)])];

    step_pts1i = interp1(step_num, step_pts1, 1:39, 'linear', 'extrap');
    step_pts2i = interp1(step_num, step_pts2, 1:39, 'linear', 'extrap');
    
    %Display the predicted step centres
    step_im = imread(['C:\isbe\density\stepwedges\images\', step_list(ii).name]);
    
    f1 = figure(...
        'windowstyle', 'normal',...
        'Units', 'centimeters ',...
        'position', [0 0 10.0 20.0],...
        'PaperPositionMode','auto');

    axes(...
        'Units', 'centimeters ',...
        'position', [0 0 10.0 20.0]);

    imagesc(step_im); axis image off; colormap(gray(256)); hold on;
    plot([step_pts1i(:,1) step_pts2i(:,1)]', [step_pts1i(:,2) step_pts2i(:,2)]');
    print('-dtiff', '-noui', '-painters', f1, '-r300', ['C:\isbe\density\stepwedges\figures\sw' zerostr(ii,2) '_steps.tif']);
    close(f1);
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%%
for ii = 1:13
    figure; imagesc(step_im); axis image; colormap(gray(256)); hold on;
    plot(x_pts, y_pts, 'rx');
    plot(x_pts(ii), y_pts(ii), 'go');
end
x_src = [39 35 31]';
y_src = [x_pts([7 9 11]), y_pts([7 9 11])];
%%
plot(x_pts(5:6), y_pts(5:6), 'r');
plot(x_pts(7:8), y_pts(7:8), 'b');
plot(x_pts(9:10), y_pts(9:10), 'g');
plot(x_pts(11:12), y_pts(11:12), 'y');
%%
s_list = dir('\\isbe-matrix\stepwedge2\stepwedge\data\*.mat');
%
wedge_vals = zeros(length(s_list), 39);
for ii = 1:length(s_list)
    load(['\\isbe-matrix\stepwedge2\stepwedge\data\' s_list(ii).name]);
    if isfield(density_data, 'wedgevals')
        wedge_vals(ii,:) = density_data.wedgevals';
    end
end
wedge_vals(~any(wedge_vals,2),:) = [];
wedge_vals(any(isnan(wedge_vals),2),:) = [];
[P_x, B_x, l_x] = princomp(wedge_vals);
%%
for ii = 1:4
    figure; 
    plot(1:39, mean(wedge_vals), '-x'); hold on;
    plot(1:39, mean(wedge_vals) + 2*sqrt(l_x(ii))*P_x(:,ii)', 'r-x');
    plot(1:39, mean(wedge_vals) - 2*sqrt(l_x(ii))*P_x(:,ii)', 'r-x');
    xlabel('Steps')
    ylabel('Step intensities');
    title(['Mean step intensities \pm 2 standard deviations of mode ' num2str(ii)]);
    saveas(gcf, ['C:\isbe\density\stepwedges\figures\mean_step_profile_mode' num2str(ii) '.bmp']);
end